#!/usr/bin/env python3
"""
Stability Monitor for RepOrtR (advanced / legacy).

This tool wraps Snakemake with additional stability features. For simple,
read-only progress dashboards, prefer `monitoring/progress_monitor.py`,
which uses the JSON files created by `workflows/progress_tracking.smk`.

Handles server disconnections and ensures assembly stability.
"""

import os
import sys
import time
import signal
import logging
import yaml
import subprocess
import psutil
from pathlib import Path
import argparse
import shutil

# Resolve configured log directory.
try:
    LOG_DIR = "logs"
    import yaml as _yaml
    with open("projects/global_config.yaml", "r") as _f:
        _cfg = _yaml.safe_load(_f) or {}
    LOG_DIR = _cfg.get("global", {}).get("log_dir", "logs")
except Exception:
    LOG_DIR = "logs"

# Setup logging
os.makedirs(LOG_DIR, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(os.path.join(LOG_DIR, 'stability_monitor.log')),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class StabilityMonitor:
    """Monitor and manage long-running assemblies"""
    
    def __init__(self, project_id, max_memory_percent=80, checkpoint_interval=3600):
        self.project_id = project_id
        self.max_memory_percent = max_memory_percent
        self.checkpoint_interval = checkpoint_interval
        self.checkpoint_file = f".snakemake/checkpoint_{project_id}.yaml"
        self.running_processes = []
        
        # Setup signal handlers
        signal.signal(signal.SIGTERM, self._signal_handler)
        signal.signal(signal.SIGINT, self._signal_handler)
        
        # Create directories
        os.makedirs(".snakemake", exist_ok=True)
        os.makedirs(LOG_DIR, exist_ok=True)
    
    def _signal_handler(self, signum, frame):
        """Handle shutdown signals gracefully"""
        logger.info(f"Received signal {signum}, initiating graceful shutdown...")
        self.save_checkpoint()
        self.cleanup()
        sys.exit(0)
    
    def save_checkpoint(self):
        """Save current state to checkpoint file"""
        checkpoint_data = {
            'timestamp': time.time(),
            'project_id': self.project_id,
            'running_processes': len(self.running_processes),
            'memory_usage': psutil.virtual_memory().percent,
            'disk_usage': psutil.disk_usage('/').percent
        }
        
        with open(self.checkpoint_file, 'w') as f:
            yaml.dump(checkpoint_data, f)
        
        logger.info(f"Checkpoint saved: {checkpoint_data}")
    
    def load_checkpoint(self):
        """Load previous checkpoint if available"""
        if os.path.exists(self.checkpoint_file):
            with open(self.checkpoint_file, 'r') as f:
                checkpoint_data = yaml.safe_load(f)
            logger.info(f"Loaded checkpoint: {checkpoint_data}")
            return checkpoint_data
        return None
    
    def monitor_memory(self):
        """Monitor memory usage and warn if high"""
        memory_percent = psutil.virtual_memory().percent
        if memory_percent > self.max_memory_percent:
            logger.warning(f"High memory usage: {memory_percent}%")
            return False
        return True
    
    def monitor_disk(self):
        """Monitor disk usage"""
        disk_percent = psutil.disk_usage('/').percent
        if disk_percent > 90:
            logger.warning(f"High disk usage: {disk_percent}%")
            return False
        return True
    
    def run_snakemake_with_stability(self, target=None, cores=16, memory="80G"):
        """Run snakemake with stability features"""
        logger.info(f"Starting Snakemake with stability features for {self.project_id}")
        
        # Load previous checkpoint
        checkpoint = self.load_checkpoint()
        
        # Start the config-driven JSON progress dashboard in the background.
        # This reads the JSON files written by `workflows/progress_tracking.smk`.
        progress_monitor_pid = None
        try:
            import subprocess
            progress_monitor_pid = subprocess.Popen(
                ['python3', '-m', 'monitoring.progress_monitor'],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            ).pid
            logger.info("Progress dashboard started")
        except Exception as e:
            logger.warning(f"Could not start progress dashboard: {e}")
        
        if checkpoint:
            logger.info("Resuming from checkpoint")
        
        # Wait a moment for beautiful progress to initialize
        time.sleep(1)
        
        # Build snakemake command with modular workflow
        cmd = [
            "snakemake",
            "-s", "Snakefile_modular",
            "--configfile", "projects/global_config.yaml",
            "--cores", str(cores),
            "--resources", f"mem_mb={int(memory.replace('G', '000'))}",
            "--rerun-incomplete",
            "--keep-going",
            "--latency-wait", "60",
            "--restart-times", "3"
        ]
        
        if target:
            # Handle multiple targets (space-separated)
            targets = target.split()
            cmd.extend(targets)
        
        # Add stability options (avoid duplicates)
        cmd.extend([
            "--keep-going",
            "--latency-wait", "60"
        ])
        
        logger.info(f"Running command: {' '.join(cmd)}")
        
        # Start monitoring in background
        self._start_monitoring()
        
        try:
            # Check for existing lock files before starting
            logger.info("🏛️ HOLY PRINCIPLE: Checking for existing lock files...")
            lock_files = [".snakemake/locks", ".snakemake/incomplete"]
            for lock_file in lock_files:
                if os.path.exists(lock_file):
                    logger.warning(f"Found existing lock file: {lock_file}")
                    logger.info("Attempting to clean up lock files...")
                    try:
                        shutil.rmtree(lock_file, ignore_errors=True)
                        logger.info(f"Cleaned up lock file: {lock_file}")
                    except Exception as e:
                        logger.warning(f"Could not clean up {lock_file}: {e}")
            
            # Run snakemake with real-time output streaming
            logger.info(f"Starting snakemake process with command: {' '.join(cmd)}")
            
            # Use Popen with real-time output streaming
            process = subprocess.Popen(
                cmd, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE,
                bufsize=1,
                universal_newlines=True
            )
            self.running_processes.append(process)
            
            logger.info(f"Snakemake process started with PID: {process.pid}")
            
            # Monitor the process with real-time output and better logging
            logger.info("🏛️ HOLY PRINCIPLE: Starting comprehensive process monitoring with real-time output...")
            monitoring_start = time.time()
            last_output_time = time.time()
            max_wait_time = 14400  # 4 hours maximum wait without output (seqclust can be very slow to start)
            
            # Set up non-blocking output reading
            import select
            import fcntl
            
            # Make stdout and stderr non-blocking
            fcntl.fcntl(process.stdout.fileno(), fcntl.F_SETFL, os.O_NONBLOCK)
            fcntl.fcntl(process.stderr.fileno(), fcntl.F_SETFL, os.O_NONBLOCK)
            
            while process.poll() is None:
                elapsed = int(time.time() - monitoring_start)
                current_time = time.time()
                
                # Check for active seqclust processes
                seqclust_running = False
                for proc in psutil.process_iter(['pid', 'cmdline']):
                    try:
                        cmdline = ' '.join(proc.info['cmdline'] or [])
                        if 'seqclust' in cmdline:
                            seqclust_running = True
                            break
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                
                # If seqclust is running, be completely hands-off and trust Snakemake's timeout
                if seqclust_running:
                    # Log progress every 120 seconds when seqclust is running
                    if elapsed % 120 == 0 and elapsed > 0:
                        logger.info(f"🏛️ HOLY PRINCIPLE: Seqclust detected - COMPLETELY HANDS-OFF monitoring. {elapsed}s elapsed, PID {process.pid} still running")
                        logger.info("🎯 Trusting Snakemake's timeout mechanism for seqclust processes - NO INTERFERENCE")
                    
                    # Don't apply ANY timeout pressure when seqclust is running
                    # Just check system resources and continue monitoring
                    if not self.monitor_memory() or not self.monitor_disk():
                        logger.warning("System resources low, saving checkpoint...")
                        self.save_checkpoint()
                    
                    time.sleep(60)  # Check every 60 seconds when seqclust is running
                    continue
                
                # Normal monitoring for non-seqclust processes
                # Log progress every 60 seconds for better responsiveness
                if elapsed % 60 == 0 and elapsed > 0:
                    logger.info(f"🏛️ HOLY PRINCIPLE: Process monitoring - {elapsed}s elapsed, PID {process.pid} still running")
                    
                    # Check for snakemake and seqclust processes
                    snakemake_count = 0
                    seqclust_count = 0
                    for proc in psutil.process_iter(['pid', 'cmdline']):
                        try:
                            cmdline = ' '.join(proc.info['cmdline'] or [])
                            if 'snakemake' in cmdline:
                                snakemake_count += 1
                            if 'seqclust' in cmdline:
                                seqclust_count += 1
                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            pass
                    
                    logger.info(f"Found {snakemake_count} active snakemake processes, {seqclust_count} seqclust processes")
                
                # Try to read output in non-blocking way
                try:
                    # Check if there's output to read
                    ready_to_read, _, _ = select.select([process.stdout, process.stderr], [], [], 0.1)
                    
                    for stream in ready_to_read:
                        line = stream.readline()
                        if line:
                            line = line.strip()
                            if line:
                                # Parse and categorize output for better logging
                                if stream == process.stdout:
                                    # Categorize stdout messages
                                    if "Building DAG" in line:
                                        logger.info(f"🏛️ HOLY PRINCIPLE: Snakemake - {line}")
                                    elif "Job" in line and "finished" in line:
                                        logger.info(f"✅ Snakemake - {line}")
                                    elif "Nothing to be done" in line:
                                        logger.info(f"ℹ️  Snakemake - {line}")
                                    elif "seqclust" in line.lower():
                                        logger.info(f"🎯 Seqclust - {line}")
                                    else:
                                        logger.info(f"Snakemake stdout: {line}")
                                else:
                                    # Categorize stderr messages
                                    if "error" in line.lower() or "failed" in line.lower():
                                        logger.error(f"❌ Snakemake error: {line}")
                                    elif "warning" in line.lower():
                                        logger.warning(f"⚠️  Snakemake warning: {line}")
                                    elif "seqclust" in line.lower():
                                        logger.info(f"🎯 Seqclust stderr: {line}")
                                    else:
                                        logger.warning(f"Snakemake stderr: {line}")
                                last_output_time = current_time
                except (OSError, IOError):
                    pass  # No output available
                
                # If no output for 300 seconds (5 minutes), log a status update (much more patient)
                if current_time - last_output_time > 300:
                    logger.info(f"🏛️ HOLY PRINCIPLE: No output for {int(current_time - last_output_time)}s, process still running...")
                    last_output_time = current_time
                
                # Timeout if no output for too long (but be very patient for seqclust)
                if current_time - last_output_time > max_wait_time:
                    logger.warning(f"🏛️ HOLY PRINCIPLE: No output for {max_wait_time}s, but continuing to monitor...")
                    # For seqclust processes, we don't kill them - they can take a very long time to start
                    # Just reset the timer and continue monitoring
                    last_output_time = current_time
                
                if not self.monitor_memory() or not self.monitor_disk():
                    logger.warning("System resources low, saving checkpoint...")
                    self.save_checkpoint()
                
                time.sleep(30)  # Check every 30 seconds for better responsiveness
            
            # Get output
            stdout, stderr = process.communicate()
            
            logger.info(f"🏛️ HOLY PRINCIPLE: Snakemake process finished with return code: {process.returncode}")
            
            # Log output with better formatting
            if stdout:
                stdout_text = stdout.strip()
                if stdout_text:
                    logger.info("🏛️ HOLY PRINCIPLE: Snakemake stdout:")
                    for line in stdout_text.split('\n'):
                        if line.strip():
                            logger.info(f"  {line.strip()}")
            
            if stderr:
                stderr_text = stderr.strip()
                if stderr_text:
                    logger.warning("🏛️ HOLY PRINCIPLE: Snakemake stderr:")
                    for line in stderr_text.split('\n'):
                        if line.strip():
                            logger.warning(f"  {line.strip()}")
            
            if process.returncode == 0:
                logger.info("🏛️ HOLY PRINCIPLE: Snakemake completed successfully")
                return True
            else:
                logger.error(f"🏛️ HOLY PRINCIPLE: Snakemake failed with return code {process.returncode}")
                
                # Provide specific guidance for common errors
                if "LockException" in stderr_text:
                    logger.error("LockException detected - this usually means another Snakemake process was running")
                    logger.error("The lock files have been cleaned up, but you may need to wait a moment before retrying")
                elif "Directory cannot be locked" in stderr_text:
                    logger.error("Directory lock issue detected - this may require manual cleanup")
                    logger.error("Try running: snakemake --unlock")
                
                return False
                
        except KeyboardInterrupt:
            logger.info("Received interrupt, saving checkpoint...")
            self.save_checkpoint()
            return False
        finally:
            self.cleanup()
    
    def _start_monitoring(self):
        """Start background monitoring"""
        def monitor_loop():
            while True:
                try:
                    self.save_checkpoint()
                    time.sleep(self.checkpoint_interval)
                except Exception as e:
                    logger.error(f"Monitoring error: {e}")
                    time.sleep(300)  # Wait 5 minutes on error
        
        import threading
        monitor_thread = threading.Thread(target=monitor_loop, daemon=True)
        monitor_thread.start()
    
    def cleanup(self):
        """Clean up resources"""
        for process in self.running_processes:
            try:
                process.terminate()
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                process.kill()
        
        self.running_processes.clear()
        
        # Stop beautiful monitoring panel if running
        try:
            subprocess.run(['pkill', '-f', 'monitoring/beautiful_progress.py'], check=False)
            logger.info("🏛️ HOLY PRINCIPLE: Beautiful monitoring panel stopped")
        except Exception as e:
            logger.warning(f"Could not stop beautiful monitoring panel: {e}")
        
        logger.info("Cleanup completed")

def main():
    parser = argparse.ArgumentParser(description='Stability Monitor for RepOrtR')
    parser.add_argument('project_id', help='Project ID to monitor')
    parser.add_argument('--target', help='Specific target to run')
    parser.add_argument('--cores', type=int, default=16, help='Number of cores')
    parser.add_argument('--memory', default='80G', help='Memory limit')
    parser.add_argument('--max-memory', type=int, default=80, help='Max memory percentage')
    parser.add_argument('--checkpoint-interval', type=int, default=3600, help='Checkpoint interval (seconds)')
    
    args = parser.parse_args()
    
    # Create monitor
    monitor = StabilityMonitor(
        args.project_id,
        max_memory_percent=args.max_memory,
        checkpoint_interval=args.checkpoint_interval
    )
    
    # Run with stability features
    success = monitor.run_snakemake_with_stability(
        target=args.target,
        cores=args.cores,
        memory=args.memory
    )
    
    if success:
        logger.info("Assembly completed successfully")
        sys.exit(0)
    else:
        logger.error("Assembly failed")
        sys.exit(1)

if __name__ == "__main__":
    main() 
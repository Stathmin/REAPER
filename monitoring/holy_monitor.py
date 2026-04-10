#!/usr/bin/env python3
"""
Holy Monitor - Comprehensive TAREAN Process Monitoring

LEGACY NOTICE:
- This script is retained for historical and experimental purposes.
- For day-to-day monitoring, prefer the config-driven dashboard in
  `monitoring/progress_monitor.py`, which reads the JSON files produced
  by `workflows/progress_tracking.smk`.

Follows holy principles for process tracking and resource monitoring.
"""

import os
import json
import time
import psutil
import glob
import subprocess
from datetime import datetime, timedelta
from pathlib import Path
import argparse


def _load_default_project_id() -> str | None:
    try:
        import yaml
        with open("projects/global_config.yaml", "r") as f:
            cfg = yaml.safe_load(f) or {}
        projects = cfg.get("projects", {}) or {}
        return next(iter(projects.keys()), None)
    except Exception:
        return None

class HolyMonitor:
    def __init__(self, project_id: str | None = None):
        self.progress_files = []
        self.seqclust_processes = {}
        self.rserv_processes = {}
        self.update_interval = 10
        self.config_file = "projects/global_config.yaml"
        self.log_dir = self._load_log_dir()
        self.project_id = project_id or _load_default_project_id() or "unknown"

    def _load_log_dir(self) -> str:
        try:
            import yaml
            with open(self.config_file, "r") as f:
                cfg = yaml.safe_load(f) or {}
            return cfg.get("global", {}).get("log_dir", "logs")
        except Exception:
            return "logs"
        
    def load_config(self):
        """Load configuration following holy principles"""
        try:
            import yaml
            with open(self.config_file, 'r') as f:
                return yaml.safe_load(f)
        except Exception as e:
            print(f"⚠️  Warning: Could not load config: {e}")
            return {}
    
    def find_progress_files(self):
        """Find all progress tracking files"""
        self.progress_files = glob.glob(f"{self.log_dir}/progress_*.json")
        
    def get_seqclust_processes(self):
        """Get all running seqclust processes with detailed info"""
        self.seqclust_processes = {}
        for proc in psutil.process_iter(['pid', 'cmdline', 'cpu_percent', 'memory_percent', 'create_time']):
            try:
                cmdline = ' '.join(proc.info['cmdline'] or [])
                if 'seqclust' in cmdline:
                    # Extract sample name from command line
                    sample = 'unknown'
                    try:
                        from post_tarean.sample_loader import load_samples_for_project
                        available_samples = list(load_samples_for_project(self.project_id).keys())
                        for sample_name in available_samples:
                            if sample_name in cmdline:
                                sample = sample_name
                                break
                    except Exception:
                        # Fallback to hardcoded check if dynamic loading fails
                        for test_sample in ['testicum', 'examicum', 'studicum']:
                            if test_sample in cmdline:
                                sample = test_sample
                                break

                    # Get detailed process info
                    proc_obj = psutil.Process(proc.info['pid'])
                    self.seqclust_processes[sample] = {
                        'pid': proc.info['pid'],
                        'cpu_percent': proc.info['cpu_percent'],
                        'memory_percent': proc.info['memory_percent'],
                        'memory_mb': proc_obj.memory_info().rss / 1024 / 1024,
                        'create_time': proc.info['create_time'],
                        'runtime': time.time() - proc.info['create_time'],
                        'cmdline': cmdline[:100] + '...' if len(cmdline) > 100 else cmdline
                    }
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
                
    def get_rserv_processes(self):
        """Get all running Rserv processes"""
        self.rserv_processes = {}
        for proc in psutil.process_iter(['pid', 'cmdline', 'cpu_percent', 'memory_percent']):
            try:
                cmdline = ' '.join(proc.info['cmdline'] or [])
                if 'Rserve' in cmdline:
                    proc_obj = psutil.Process(proc.info['pid'])
                    self.rserv_processes[proc.info['pid']] = {
                        'pid': proc.info['pid'],
                        'cpu_percent': proc.info['cpu_percent'],
                        'memory_percent': proc.info['memory_percent'],
                        'memory_mb': proc_obj.memory_info().rss / 1024 / 1024,
                        'cmdline': cmdline[:100] + '...' if len(cmdline) > 100 else cmdline
                    }
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
                
    def load_progress_data(self, progress_file):
        """Load progress data from JSON file"""
        try:
            with open(progress_file, 'r') as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return None
            
    def format_progress_bar(self, progress, width=20):
        """Format a beautiful progress bar"""
        filled = int(width * progress / 100)
        bar = '█' * filled + '░' * (width - filled)
        return f"[{bar}] {progress:.1f}%"
        
    def format_time(self, seconds):
        """Format time in human-readable format"""
        if seconds is None:
            return "Unknown"
        minutes = int(seconds // 60)
        secs = int(seconds % 60)
        if minutes > 60:
            hours = minutes // 60
            minutes = minutes % 60
            return f"{hours}h {minutes}m {secs}s"
        else:
            return f"{minutes}m {secs}s"
            
    def calculate_eta(self, progress_data):
        """Calculate ETA based on holy principles"""
        if not progress_data:
            return None
            
        config = self.load_config()
        stages_config = config.get('global', {}).get('progress_tracking', {}).get('stages', {})
        
        current_stage = progress_data.get('current_stage', 'initialization')
        current_stage_data = progress_data['progress'].get(current_stage, {})
        
        if current_stage_data.get('status') != 'running':
            return None
            
        # Calculate remaining time for current stage
        if current_stage_data.get('start_time'):
            elapsed = time.time() - current_stage_data['start_time']
            stage_duration = stages_config.get(current_stage, {}).get('estimated_duration', 300)
            remaining_current = max(0, stage_duration - elapsed)
        else:
            remaining_current = stages_config.get(current_stage, {}).get('estimated_duration', 300)
        
        # Add time for remaining stages
        total_remaining = remaining_current
        for stage_name, stage_data in progress_data['progress'].items():
            if stage_data['status'] == 'pending':
                stage_duration = stages_config.get(stage_name, {}).get('estimated_duration', 300)
                total_remaining += stage_duration
                
        return total_remaining
            
    def display_sample_progress(self, sample, progress_data):
        """Display progress for a single sample"""
        print(f"\n🎯 TAREAN Progress: {sample.upper()}")
        print("=" * 60)
        
        # Overall progress
        overall_progress = 0
        completed_stages = 0
        total_stages = len(progress_data['progress'])
        
        for stage_name, stage_data in progress_data['progress'].items():
            if stage_data['status'] == 'completed':
                completed_stages += 1
                overall_progress += 100
            elif stage_data['status'] == 'running':
                overall_progress += stage_data['progress']
                
        overall_progress /= total_stages
        
        print(f"📊 Overall Progress: {overall_progress:.1f}%")
        print(f"⏱️  Current Stage: {progress_data['current_stage'].replace('_', ' ').title()}")
        
        # Resource usage
        if sample in self.seqclust_processes:
            proc_info = self.seqclust_processes[sample]
            print(f"💻 CPU Usage: {proc_info['cpu_percent']:.1f}%")
            print(f"💾 Memory Usage: {proc_info['memory_mb']:.1f} MB")
            print(f"⏱️  Runtime: {self.format_time(proc_info['runtime'])}")
        else:
            print("💻 CPU Usage: N/A (process not found)")
            print("💾 Memory Usage: N/A")
            
        print(f"\n📋 Stage Status:")
        
        # Display each stage
        for stage_name, stage_data in progress_data['progress'].items():
            status_icon = '✅' if stage_data['status'] == 'completed' else '🔄' if stage_data['status'] == 'running' else '⏳'
            stage_display = stage_name.replace('_', ' ').title()
            progress_bar = self.format_progress_bar(stage_data['progress'])
            print(f"  {status_icon} {stage_display}: {progress_bar}")
            
        # ETA
        eta = self.calculate_eta(progress_data)
        if eta:
            eta_time = self.format_time(eta)
            print(f"\n⏰ Estimated Time Remaining: {eta_time}")
            
        print(f"\n🕐 Last Updated: {datetime.now().strftime('%H:%M:%S')}")
        
    def display_process_balance(self):
        """Display process balance and resource distribution"""
        print(f"\n⚖️  Process Balance Analysis")
        print("=" * 60)
        
        if not self.seqclust_processes:
            print("⏸️  No seqclust processes running")
            return
            
        total_cpu = sum(proc['cpu_percent'] for proc in self.seqclust_processes.values())
        total_memory = sum(proc['memory_mb'] for proc in self.seqclust_processes.values())
        
        print(f"📊 Total CPU Usage: {total_cpu:.1f}%")
        print(f"💾 Total Memory Usage: {total_memory:.1f} MB")
        print(f"🔄 Active Processes: {len(self.seqclust_processes)}")
        
        # Check if processes are balanced
        if len(self.seqclust_processes) > 1:
            cpu_values = [proc['cpu_percent'] for proc in self.seqclust_processes.values()]
            cpu_variance = max(cpu_values) - min(cpu_values)
            
            if cpu_variance < 10:
                print("✅ CPU usage is well balanced")
            elif cpu_variance < 30:
                print("⚠️  CPU usage is moderately balanced")
            else:
                print("❌ CPU usage is poorly balanced")
                
        # Display individual process info
        for sample, proc_info in self.seqclust_processes.items():
            print(f"  📊 {sample}: CPU {proc_info['cpu_percent']:.1f}%, MEM {proc_info['memory_mb']:.1f}MB")
            
    def display_output_tracking(self):
        """Display output file tracking"""
        print(f"\n📁 Output Tracking")
        print("=" * 60)
        
        try:
            from post_tarean.sample_loader import load_samples_for_project
            samples = list(load_samples_for_project(self.project_id).keys())
        except Exception as e:
            print(f"Warning: Could not load samples dynamically: {e}")
            samples = []
        
        for sample in samples:
            tarean_dir = f"projects/{self.project_id}/samples/{sample}/tarean"
            done_file = f"{tarean_dir}/tarean.done"
            log_file = f"{self.log_dir}/seqclust_{self.project_id}_{sample}.log"
            
            status_icon = '✅' if os.path.exists(done_file) else '⏳' if os.path.exists(log_file) else '❌'
            print(f"  {status_icon} {sample}:")
            print(f"    📄 Done file: {'✅' if os.path.exists(done_file) else '❌'}")
            print(f"    📝 Log file: {'✅' if os.path.exists(log_file) else '❌'}")
            
            if os.path.exists(log_file):
                log_size = os.path.getsize(log_file)
                print(f"    📏 Log size: {log_size} bytes")
                
    def display_dashboard(self):
        """Display the main holy monitoring dashboard"""
        os.system('clear')
        print("🏛️  Holy RepOrtR Monitor")
        print("=" * 80)
        print(f"📅 {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 80)
        
        # Update data
        self.find_progress_files()
        self.get_seqclust_processes()
        self.get_rserv_processes()
        
        # Display progress for each sample
        for progress_file in self.progress_files:
            progress_data = self.load_progress_data(progress_file)
            if progress_data:
                sample = progress_data['sample']
                self.display_sample_progress(sample, progress_data)
                
        # Display process balance
        self.display_process_balance()
        
        # Display output tracking
        self.display_output_tracking()
        
        # System resources
        print(f"\n💻 System Resources:")
        cpu_percent = psutil.cpu_percent(interval=1)
        memory = psutil.virtual_memory()
        print(f"  CPU Usage: {cpu_percent:.1f}%")
        print(f"  Memory Usage: {memory.percent:.1f}% ({memory.used // (1024**3):.1f}GB / {memory.total // (1024**3):.1f}GB)")
        
        # Rserv processes
        if self.rserv_processes:
            print(f"\n🖥️  Rserv Processes ({len(self.rserv_processes)}):")
            for pid, proc_info in self.rserv_processes.items():
                print(f"  PID {pid}: CPU {proc_info['cpu_percent']:.1f}%, MEM {proc_info['memory_mb']:.1f}MB")
        
        print(f"\n🔄 Auto-refresh every {self.update_interval} seconds (Press Ctrl+C to exit)")
        
    def run(self):
        """Run the holy monitor"""
        try:
            while True:
                self.display_dashboard()
                time.sleep(self.update_interval)
        except KeyboardInterrupt:
            print("\n👋 Holy monitoring stopped")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Legacy holy monitor")
    parser.add_argument("--project-id", help="Project ID (defaults to first in global_config.yaml)")
    args = parser.parse_args()
    monitor = HolyMonitor(project_id=args.project_id)
    monitor.run()

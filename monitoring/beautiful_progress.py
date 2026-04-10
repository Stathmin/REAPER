#!/usr/bin/env python3
"""
Beautiful Progress Display for RepOrtR
Real-time seqclust progress tracking with beautiful output
"""

import os
import json
import time
import psutil
import glob
from datetime import datetime, timedelta
from pathlib import Path
import argparse


def _load_default_project_id():
    try:
        import yaml

        with open("projects/global_config.yaml", "r") as f:
            cfg = yaml.safe_load(f) or {}
        projects = cfg.get("projects", {}) or {}
        return next(iter(projects.keys()), None)
    except Exception:
        return None

class BeautifulProgress:
    def __init__(self, project_id=None):
        self.progress_files = []
        self.seqclust_processes = {}
        self.update_interval = 10
        self.log_dir = self._load_log_dir()
        self.project_id = project_id or _load_default_project_id() or "unknown"
        
    def _load_log_dir(self) -> str:
        try:
            import yaml
            with open("projects/global_config.yaml", "r") as f:
                cfg = yaml.safe_load(f) or {}
            return cfg.get("global", {}).get("log_dir", "logs")
        except Exception:
            return "logs"
        
    def find_progress_files(self):
        """Find all progress tracking files"""
        self.progress_files = glob.glob(f"{self.log_dir}/progress_*.json")
        
    def get_seqclust_processes(self):
        """Get all running seqclust processes"""
        self.seqclust_processes = {}
        for proc in psutil.process_iter(['pid', 'cmdline', 'cpu_percent', 'memory_percent']):
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
                        if 'examicum' in cmdline:
                            sample = 'examicum'
                        elif 'studicum' in cmdline:
                            sample = 'studicum'
                        elif 'testicum' in cmdline:
                            sample = 'testicum'
                    
                    self.seqclust_processes[sample] = {
                        'pid': proc.info['pid'],
                        'cpu_percent': proc.info['cpu_percent'],
                        'memory_percent': proc.info['memory_percent'],
                        'memory_mb': proc.memory_info().rss / 1024 / 1024
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
        if seconds is None or seconds <= 0:
            return "Unknown"
        
        # Use actual seconds, but cap display at 24 hours for readability
        display_seconds = min(seconds, 86400)
        
        minutes = int(display_seconds // 60)
        secs = int(display_seconds % 60)
        
        if minutes >= 1440:  # 24 hours
            days = minutes // 1440
            hours = (minutes % 1440) // 60
            minutes = minutes % 60
            display_str = f"{days}d {hours}h {minutes}m"
            if seconds > 86400:
                return f"{display_str} (actual: {seconds//3600}h)"
            return display_str
        elif minutes >= 60:
            hours = minutes // 60
            minutes = minutes % 60
            return f"{hours}h {minutes}m {secs}s"
        else:
            return f"{minutes}m {secs}s"
            
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
        if progress_data.get('eta'):
            eta_time = self.format_time(progress_data['eta'])
            print(f"\n⏰ Estimated Time Remaining: {eta_time}")
            
        # Show recent seqclust output summary
        project = progress_data.get("project") or self.project_id
        self.display_seqclust_summary(project, sample)
            
        print(f"\n🕐 Last Updated: {datetime.now().strftime('%H:%M:%S')}")
    
    def display_seqclust_summary(self, project, sample):
        """Display recent seqclust output summary"""
        try:
            log_file = f"{self.log_dir}/seqclust_{project}_{sample}.log"
            if os.path.exists(log_file):
                # Get last few lines of seqclust output
                with open(log_file, 'r') as f:
                    lines = f.readlines()
                    if lines:
                        # Show last 3 relevant lines
                        recent_lines = []
                        for line in reversed(lines[-10:]):  # Check last 10 lines
                            line = line.strip()
                            if line and any(keyword in line.lower() for keyword in ['clustering', 'assembly', 'annotation', 'progress', 'stage']):
                                recent_lines.append(line)
                                if len(recent_lines) >= 3:
                                    break
                        
                        if recent_lines:
                            print(f"\n🎯 Recent Seqclust Activity:")
                            for line in recent_lines:
                                # Truncate long lines for display
                                if len(line) > 80:
                                    line = line[:77] + "..."
                                print(f"  {line}")
        except Exception as e:
            pass  # Ignore errors in summary display
        
    def display_dashboard(self):
        """Display the main progress dashboard"""
        os.system('clear')
        print("🎯 RepOrtR Progress Dashboard")
        print("=" * 80)
        print(f"📅 {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 80)
        
        # Update data
        self.find_progress_files()
        self.get_seqclust_processes()
        
        if not self.progress_files:
            print("\n🏛️ HOLY PRINCIPLE: RepOrtR Progress Dashboard")
            print("=" * 80)
            print(f"📅 {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            print("=" * 80)
            
            # Check if any seqclust processes are running
            if self.seqclust_processes:
                print(f"\n🔄 TAREAN Analysis in Progress:")
                for sample, proc_info in self.seqclust_processes.items():
                    print(f"  🎯 {sample.upper()}: PID {proc_info['pid']}")
                    print(f"     CPU: {proc_info['cpu_percent']:.1f}% | Memory: {proc_info['memory_mb']:.1f} MB")
                print(f"\n⏳ Progress tracking files will be created as TAREAN analysis progresses...")
            else:
                print(f"\n⏳ Waiting for TAREAN analysis to begin...")
                print(f"   Progress tracking will start when seqclust processes are detected")
                print(f"   Checking for snakemake processes that may be building dependencies...")
                
                # Check for snakemake processes
                snakemake_processes = []
                for proc in psutil.process_iter(['pid', 'cmdline', 'cpu_percent']):
                    try:
                        cmdline = ' '.join(proc.info['cmdline'] or [])
                        if 'snakemake' in cmdline:
                            snakemake_processes.append({
                                'pid': proc.info['pid'],
                                'cpu_percent': proc.info['cpu_percent']
                            })
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                
                if snakemake_processes:
                    print(f"   Found {len(snakemake_processes)} snakemake processes:")
                    for proc in snakemake_processes:
                        print(f"     PID {proc['pid']}: CPU {proc['cpu_percent']:.1f}%")
                    print(f"   These may be building dependencies before seqclust starts...")
            
            # System resources
            print(f"\n💻 System Resources:")
            cpu_percent = psutil.cpu_percent(interval=0)  # Non-blocking call
            memory = psutil.virtual_memory()
            print(f"  CPU Usage: {cpu_percent:.1f}%")
            print(f"  Memory Usage: {memory.percent:.1f}% ({memory.used // (1024**3):.1f}GB / {memory.total // (1024**3):.1f}GB)")
            
            print(f"\n🔄 Auto-refresh every {self.update_interval} seconds (Press Ctrl+C to exit)")
            return
            
        # Display progress for each sample
        for progress_file in self.progress_files:
            progress_data = self.load_progress_data(progress_file)
            if progress_data:
                sample = progress_data['sample']
                self.display_sample_progress(sample, progress_data)
                
        # System resources
        print(f"\n💻 System Resources:")
        cpu_percent = psutil.cpu_percent(interval=0)  # Non-blocking call
        memory = psutil.virtual_memory()
        print(f"  CPU Usage: {cpu_percent:.1f}%")
        print(f"  Memory Usage: {memory.percent:.1f}% ({memory.used // (1024**3):.1f}GB / {memory.total // (1024**3):.1f}GB)")
        
        print(f"\n🔄 Auto-refresh every {self.update_interval} seconds (Press Ctrl+C to exit)")
        
    def run(self):
        """Run the beautiful progress display"""
        try:
            while True:
                self.display_dashboard()
                time.sleep(self.update_interval)
        except KeyboardInterrupt:
            print("\n👋 Progress monitoring stopped")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Legacy progress dashboard (process-scanning)")
    parser.add_argument("--project-id", help="Project ID (defaults to first in global_config.yaml)")
    args = parser.parse_args()
    progress = BeautifulProgress(project_id=args.project_id)
    progress.run()

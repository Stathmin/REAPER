#!/usr/bin/env python3
"""
Pipeline Monitoring Script (legacy).

LEGACY NOTICE:
- This script predates the JSON-based progress tracking rules in
  `workflows/progress_tracking.smk`.
- For current projects, use `monitoring/progress_monitor.py`, which is
  driven by `projects/global_config.yaml` and the standard progress files.

Monitors TAREAN processes and resource usage.
"""

import os
import time
import subprocess
import psutil
from datetime import datetime
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

def get_seqclust_processes():
    """Get all running seqclust processes"""
    processes = []
    for proc in psutil.process_iter(['pid', 'name', 'cmdline', 'cpu_percent', 'memory_percent']):
        try:
            if 'seqclust' in ' '.join(proc.info['cmdline'] or []):
                processes.append(proc.info)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return processes

def get_rserv_processes():
    """Get all running Rserv processes"""
    processes = []
    for proc in psutil.process_iter(['pid', 'name', 'cmdline', 'cpu_percent', 'memory_percent']):
        try:
            if 'Rserve' in ' '.join(proc.info['cmdline'] or []):
                processes.append(proc.info)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return processes

def check_completion_files(project_id: str):
    """Check completion status of samples"""
    try:
        from post_tarean.sample_loader import load_samples_for_project
        samples = list(load_samples_for_project(project_id).keys())
    except Exception as e:
        print(f"Warning: Could not load samples dynamically: {e}")
        samples = []  # Fallback to empty list
    
    status = {}
    
    for sample in samples:
        done_file = f"projects/{project_id}/samples/{sample}/tarean/tarean.done"
        status[sample] = os.path.exists(done_file)
    
    return status

def main(project_id: str):
    """Main monitoring function"""
    print("🔍 RepOrtR Pipeline Monitor")
    print("=" * 50)
    
    while True:
        # Clear screen
        os.system('clear')
        
        print(f"📊 Pipeline Status - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 50)
        
        # Check completion status
        print("\n📋 Sample Completion Status:")
        status = check_completion_files(project_id)
        for sample, completed in status.items():
            status_icon = "✅" if completed else "⏳"
            print(f"  {status_icon} {sample}: {'Completed' if completed else 'Running/Pending'}")
        
        # Check running processes
        print("\n🔄 Active Processes:")
        seqclust_procs = get_seqclust_processes()
        rserv_procs = get_rserv_processes()
        
        if seqclust_procs:
            print("  📊 TAREAN/seqclust processes:")
            for proc in seqclust_procs:
                cmd = ' '.join(proc['cmdline'] or [])
                # Try to extract sample name from command
                sample = 'unknown'
                try:
                    from post_tarean.sample_loader import load_samples_for_project
                    available_samples = list(load_samples_for_project(project_id).keys())
                    for sample_name in available_samples:
                        if sample_name in cmd:
                            sample = sample_name
                            break
                except Exception:
                    # Fallback to hardcoded check if dynamic loading fails
                    if 'examicum' in cmd:
                        sample = 'examicum'
                    elif 'studicum' in cmd:
                        sample = 'studicum'
                print(f"    PID {proc['pid']}: {sample} (CPU: {proc['cpu_percent']:.1f}%, MEM: {proc['memory_percent']:.1f}%)")
        else:
            print("  ⏸️  No TAREAN processes running")
        
        if rserv_procs:
            print("  🖥️  Rserv processes:")
            for proc in rserv_procs:
                print(f"    PID {proc['pid']}: Rserve (CPU: {proc['cpu_percent']:.1f}%, MEM: {proc['memory_percent']:.1f}%)")
        else:
            print("  ⏸️  No Rserv processes running")
        
        # System resources
        print("\n💻 System Resources:")
        cpu_percent = psutil.cpu_percent(interval=1)
        memory = psutil.virtual_memory()
        print(f"  CPU Usage: {cpu_percent:.1f}%")
        print(f"  Memory Usage: {memory.percent:.1f}% ({memory.used // (1024**3):.1f}GB / {memory.total // (1024**3):.1f}GB)")
        
        print("\n" + "=" * 50)
        print("Press Ctrl+C to exit")
        
        time.sleep(5)

if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(description="Legacy pipeline monitor")
        parser.add_argument("--project-id", help="Project ID (defaults to first in global_config.yaml)")
        args = parser.parse_args()
        pid = args.project_id or _load_default_project_id() or "unknown"
        main(pid)
    except KeyboardInterrupt:
        print("\n👋 Monitoring stopped")

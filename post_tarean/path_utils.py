#!/usr/bin/env python3
"""
Path Utilities for RepOrtR

This module provides utilities for handling test vs production paths properly.

Author: RepOrtR Team
Date: 2025
"""

import os
from pathlib import Path
from typing import Union

def is_test_project(project_id: str) -> bool:
    """
    Check if a project ID is a test project
    
    Args:
        project_id: Project identifier
        
    Returns:
        True if test project, False otherwise
    """
    test_indicators = ['toy_project', 'test_', 'mock_', 'dummy_', 'example_']
    return any(indicator in project_id.lower() for indicator in test_indicators)

def get_test_output_dir(tool_name: str) -> Path:
    """
    Get test output directory for a specific tool
    
    Args:
        tool_name: Name of the tool (e.g., 'ncbi_gathering', 'blast_comparison')
        
    Returns:
        Path to test output directory
    """
    return Path(f"test_outputs/{tool_name}")

def get_project_path(project_id: str, tool_name: str = None) -> Path:
    """
    Get appropriate project path (test or production)
    
    Args:
        project_id: Project identifier
        tool_name: Tool name for test outputs
        
    Returns:
        Path to project directory
    """
    if is_test_project(project_id):
        if tool_name:
            return get_test_output_dir(tool_name)
        else:
            return Path("test_outputs")
    else:
        return Path(f"projects/{project_id}")

def ensure_output_dir(path: Union[str, Path], create: bool = True) -> Path:
    """
    Ensure output directory exists
    
    Args:
        path: Directory path
        create: Whether to create directory if it doesn't exist
        
    Returns:
        Path object
    """
    path_obj = Path(path)
    if create:
        path_obj.mkdir(parents=True, exist_ok=True)
    return path_obj

def get_output_path(project_id: str, tool_name: str, filename: str) -> Path:
    """
    Get output path for a file, respecting test vs production projects
    
    Args:
        project_id: Project identifier
        tool_name: Tool name
        filename: Filename
        
    Returns:
        Path to output file
    """
    if is_test_project(project_id):
        output_dir = get_test_output_dir(tool_name)
    else:
        output_dir = Path(f"projects/{project_id}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir / filename

def cleanup_test_outputs():
    """Clean up test output directories"""
    test_outputs = Path("test_outputs")
    if test_outputs.exists():
        import shutil
        shutil.rmtree(test_outputs)
        print("🧹 Cleaned up test outputs")

def list_test_outputs():
    """List all test output directories and their contents"""
    test_outputs = Path("test_outputs")
    if not test_outputs.exists():
        print("No test outputs directory found")
        return
    
    print("📁 Test Outputs Directory Structure:")
    print("=" * 50)
    
    for tool_dir in test_outputs.iterdir():
        if tool_dir.is_dir():
            print(f"\n🔧 {tool_dir.name}/")
            for item in tool_dir.iterdir():
                if item.is_file():
                    size = item.stat().st_size
                    print(f"  📄 {item.name} ({size} bytes)")
                elif item.is_dir():
                    file_count = len(list(item.glob("*")))
                    print(f"  📁 {item.name}/ ({file_count} items)")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        command = sys.argv[1]
        
        if command == "cleanup":
            cleanup_test_outputs()
        elif command == "list":
            list_test_outputs()
        else:
            print("Usage: python path_utils.py [cleanup|list]")
    else:
        print("Path Utilities for RepOrtR")
        print("Available commands: cleanup, list")


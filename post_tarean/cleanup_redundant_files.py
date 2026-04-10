#!/usr/bin/env python3
"""
Cleanup script for redundant files after migration to consolidated pipeline
This script helps safely remove the old redundant files
"""

import os
import shutil
import argparse
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class RedundantFileCleaner:
    """Cleaner for redundant files after pipeline consolidation"""
    
    def __init__(self, backup_dir: str = None):
        self.backup_dir = backup_dir or "backup_redundant_files"
        self.redundant_files = {
            # BLAST variants (consolidated into blast_consolidated.py)
            'blast_variants': [
                'blast.py',
                'blast_safe.py',
                'blast_Kroupin.py',
                'blastold.py',
                'super_blast.py'
            ],
            
            # Representation variants (consolidated into pipeline.py)
            'representation_variants': [
                'represenation.py',
                'represenation_safe.py',
                'represenation_Kroupin.py',
                'represenation_Kroupin_abomination.py',
                'represenation (copy).py'
            ],
            
            # Word processing (consolidated into pipeline.py)
            'word_processing': [
                'worder.py',
                'worder_new.py'
            ],
            
            # Analysis scripts (consolidated into pipeline.py)
            'analysis_scripts': [
                'annotate_repeatome_structure.py',
                'repeat_sums_by_group.py',
                'types.py'
            ],
            
            # Utility scripts (consolidated into utils.py)
            'utility_scripts': [
                'rm_homology.py',
                'screenshot_repex.py',
                'katya_csv.py',
                'kate_rep.py'
            ]
        }
        
        # Core files that should NOT be removed
        self.core_files = [
            'pipeline.py',
            'blast_consolidated.py',
            'utils.py',
            'config.py',
            'README.md',
            'cleanup_redundant_files.py'
        ]
    
    def list_redundant_files(self, dry_run: bool = True) -> dict:
        """List all redundant files found in the directory"""
        found_files = {}
        current_dir = Path('.')
        
        for category, files in self.redundant_files.items():
            found_files[category] = []
            for file in files:
                file_path = current_dir / file
                if file_path.exists():
                    found_files[category].append(str(file_path))
        
        # Report findings
        total_found = sum(len(files) for files in found_files.values())
        logger.info(f"Found {total_found} redundant files:")
        
        for category, files in found_files.items():
            if files:
                logger.info(f"  {category}: {len(files)} files")
                for file in files:
                    action = "Would remove" if dry_run else "Removing"
                    logger.info(f"    {action}: {file}")
        
        return found_files
    
    def create_backup(self, files_to_backup: list) -> bool:
        """Create backup of files before removal"""
        if not files_to_backup:
            return True
        
        try:
            # Create backup directory
            backup_path = Path(self.backup_dir)
            backup_path.mkdir(exist_ok=True)
            
            logger.info(f"Creating backup in {backup_path}")
            
            # Copy files to backup
            for file_path in files_to_backup:
                src = Path(file_path)
                dst = backup_path / src.name
                
                if src.exists():
                    shutil.copy2(src, dst)
                    logger.info(f"  Backed up: {file_path}")
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to create backup: {e}")
            return False
    
    def remove_files(self, files_to_remove: list, force: bool = False) -> bool:
        """Remove redundant files"""
        if not files_to_remove:
            return True
        
        failed_removals = []
        
        for file_path in files_to_remove:
            try:
                path = Path(file_path)
                if path.exists():
                    if force or self._confirm_removal(str(path)):
                        path.unlink()
                        logger.info(f"Removed: {file_path}")
                    else:
                        logger.info(f"Skipped: {file_path}")
                else:
                    logger.warning(f"File not found: {file_path}")
                    
            except Exception as e:
                logger.error(f"Failed to remove {file_path}: {e}")
                failed_removals.append(file_path)
        
        if failed_removals:
            logger.error(f"Failed to remove {len(failed_removals)} files")
            return False
        
        return True
    
    def _confirm_removal(self, file_path: str) -> bool:
        """Confirm file removal (for interactive mode)"""
        response = input(f"Remove {file_path}? (y/N): ").strip().lower()
        return response in ['y', 'yes']
    
    def cleanup(self, dry_run: bool = True, backup: bool = True, 
                force: bool = False, interactive: bool = False) -> bool:
        """Main cleanup function"""
        logger.info("Starting redundant file cleanup...")
        
        # List redundant files
        found_files = self.list_redundant_files(dry_run)
        all_files = []
        for files in found_files.values():
            all_files.extend(files)
        
        if not all_files:
            logger.info("No redundant files found to clean up")
            return True
        
        if dry_run:
            logger.info("DRY RUN MODE - No files will be modified")
            return True
        
        # Create backup if requested
        if backup and all_files:
            if not self.create_backup(all_files):
                logger.error("Backup failed, aborting cleanup")
                return False
        
        # Remove files
        if interactive:
            return self.remove_files(all_files, force=False)
        else:
            return self.remove_files(all_files, force=force)
    
    def verify_cleanup(self) -> dict:
        """Verify that cleanup was successful"""
        logger.info("Verifying cleanup...")
        
        verification = {
            'core_files_present': [],
            'core_files_missing': [],
            'redundant_files_remaining': [],
            'backup_created': False
        }
        
        # Check core files
        for file in self.core_files:
            if Path(file).exists():
                verification['core_files_present'].append(file)
            else:
                verification['core_files_missing'].append(file)
        
        # Check for remaining redundant files
        for category, files in self.redundant_files.items():
            for file in files:
                if Path(file).exists():
                    verification['redundant_files_remaining'].append(file)
        
        # Check backup
        if Path(self.backup_dir).exists():
            verification['backup_created'] = True
        
        # Report results
        logger.info(f"Core files present: {len(verification['core_files_present'])}")
        logger.info(f"Core files missing: {len(verification['core_files_missing'])}")
        logger.info(f"Redundant files remaining: {len(verification['redundant_files_remaining'])}")
        logger.info(f"Backup created: {verification['backup_created']}")
        
        return verification
    
    def restore_from_backup(self) -> bool:
        """Restore files from backup"""
        backup_path = Path(self.backup_dir)
        
        if not backup_path.exists():
            logger.error(f"Backup directory not found: {self.backup_dir}")
            return False
        
        logger.info(f"Restoring files from {backup_path}")
        
        restored_count = 0
        for backup_file in backup_path.iterdir():
            if backup_file.is_file():
                try:
                    shutil.copy2(backup_file, backup_file.name)
                    logger.info(f"Restored: {backup_file.name}")
                    restored_count += 1
                except Exception as e:
                    logger.error(f"Failed to restore {backup_file.name}: {e}")
        
        logger.info(f"Restored {restored_count} files")
        return True


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Clean up redundant files after migrating to consolidated pipeline'
    )
    parser.add_argument(
        '--dry-run', 
        action='store_true',
        help='Show what would be removed without actually removing files'
    )
    parser.add_argument(
        '--no-backup',
        action='store_true',
        help='Skip creating backup before removal'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Force removal without confirmation'
    )
    parser.add_argument(
        '--interactive',
        action='store_true',
        help='Interactive mode - confirm each file removal'
    )
    parser.add_argument(
        '--verify',
        action='store_true',
        help='Verify cleanup results'
    )
    parser.add_argument(
        '--restore',
        action='store_true',
        help='Restore files from backup'
    )
    parser.add_argument(
        '--backup-dir',
        default='backup_redundant_files',
        help='Backup directory name (default: backup_redundant_files)'
    )
    
    args = parser.parse_args()
    
    # Initialize cleaner
    cleaner = RedundantFileCleaner(args.backup_dir)
    
    if args.restore:
        success = cleaner.restore_from_backup()
        exit(0 if success else 1)
    
    if args.verify:
        verification = cleaner.verify_cleanup()
        exit(0 if not verification['redundant_files_remaining'] else 1)
    
    # Perform cleanup
    success = cleaner.cleanup(
        dry_run=args.dry_run,
        backup=not args.no_backup,
        force=args.force,
        interactive=args.interactive
    )
    
    if success and not args.dry_run:
        logger.info("Cleanup completed successfully")
        cleaner.verify_cleanup()
    
    exit(0 if success else 1)


if __name__ == "__main__":
    main() 
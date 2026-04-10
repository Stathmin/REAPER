#!/usr/bin/env python3
"""
REAPER Installation Module
Creates conda environments from repository YAMLs:
- reportr: Snakemake + Python tooling used by the workflow and post-processing
- repeatexplorer: RepeatExplorer/TAREAN runtime for building and running `repex_tarean/seqclust`

This installer no longer creates a repo-root `seqclust` wrapper. Snakemake is expected
to run rules under the appropriate `conda:` environments.
"""

import os
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path

# Consistent formatting constants
SUCCESS = "✅"
ERROR = "❌"
WARNING = "⚠️"
INFO = "ℹ️"

def print_header(title):
    """Print consistent header"""
    print(f"\n{'='*60}")
    print(f"{title}")
    print(f"{'='*60}")

def print_section(title):
    """Print consistent section header"""
    print(f"\n--- {title} ---")

def run_command(cmd, description, check=True, cwd=None):
    """Run a command with consistent formatting"""
    print(f"Running: {description}")
    print(f"Command: {' '.join(cmd)}")
    if cwd:
        print(f"Working directory: {cwd}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, cwd=cwd)
        
        if result.returncode == 0:
            print(f"{SUCCESS} {description} completed successfully")
            return True
        else:
            print(f"{ERROR} {description} failed")
            if result.stderr:
                print(f"Error: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"{ERROR} {description} timed out")
        return False
    except Exception as e:
        print(f"{ERROR} {description} error: {e}")
        return False

def find_conda():
    """Find conda installation in common locations"""
    print_section("Finding Conda Installation")
    
    # Common conda locations
    conda_locations = [
        "conda",
        "~/miniconda3/bin/conda",
        "~/anaconda3/bin/conda",
        "~/miniconda/bin/conda",
        "~/anaconda/bin/conda",
        "/opt/conda/bin/conda",
        "/usr/local/bin/conda",
        "/opt/miniconda3/bin/conda",
        "/opt/anaconda3/bin/conda"
    ]
    
    # Expand home directory
    home = os.path.expanduser("~")
    conda_locations = [loc.replace("~", home) for loc in conda_locations]
    
    for conda_path in conda_locations:
        try:
            result = subprocess.run([conda_path, "--version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                print(f"{SUCCESS} Found conda at: {conda_path}")
                print(f"   Version: {result.stdout.strip()}")
                return conda_path
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    
    return None

def check_conda_available():
    """Check if conda is available and provide setup guidance"""
    print_section("Checking Conda Availability")
    
    # First try PATH
    try:
        result = subprocess.run(["conda", "--version"], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            print(f"{SUCCESS} Conda is available in PATH: {result.stdout.strip()}")
            return "conda"
    except FileNotFoundError:
        pass
    
    # Try to find conda in common locations
    conda_path = find_conda()
    if conda_path:
        print(f"{INFO} Using conda from: {conda_path}")
        return conda_path
    
    # Conda not found - provide installation guidance
    print(f"{ERROR} Conda not found in PATH or common locations")
    print(f"\n{INFO} Conda Installation Options:")
    print(f"  1. Install Miniconda (recommended):")
    print(f"     - Download from: https://docs.conda.io/en/latest/miniconda.html")
    print(f"     - Run installer and restart terminal")
    print(f"     - Or run: source ~/.bashrc (or ~/.zshrc)")
    print(f"  2. Install Anaconda:")
    print(f"     - Download from: https://www.anaconda.com/products/distribution")
    print(f"     - Run installer and restart terminal")
    print(f"  3. If conda is already installed but not in PATH:")
    print(f"     - Add to PATH: export PATH=\"~/miniconda3/bin:$PATH\"")
    print(f"     - Or initialize: conda init")
    print(f"     - Then restart terminal or run: source ~/.bashrc")
    
    print(f"\n{INFO} After installing conda, run this script again.")
    return None

def _ensure_env_from_yaml(conda_path: str, env_name: str, env_yaml: str) -> bool:
    """Create or update a conda env from a YAML file."""
    print_section(f"Ensuring conda env: {env_name}")

    if not os.path.exists(env_yaml):
        print(f"{ERROR} Missing environment YAML: {env_yaml}")
        return False

    # Check if environment exists
    env_exists = False
    try:
        result = subprocess.run([conda_path, "env", "list"], capture_output=True, text=True, timeout=30)
        if result.returncode == 0 and env_name in result.stdout:
            env_exists = True
    except Exception as e:
        print(f"{WARNING} Could not check existing environments: {e}")

    if env_exists:
        # Prefer update-in-place; fall back to recreate if update fails.
        if run_command([conda_path, "env", "update", "-n", env_name, "-f", env_yaml], f"Updating {env_name} from {env_yaml}"):
            return True
        print(f"{WARNING} Update failed; will try recreate for {env_name}")
        response = input(f"Remove and recreate '{env_name}'? (y/N): ")
        if response.lower() not in ["y", "yes"]:
            print(f"{INFO} Using existing environment")
            return True
        if not run_command([conda_path, "env", "remove", "-n", env_name], f"Removing {env_name}", check=False):
            print(f"{WARNING} Normal removal failed; you may need to remove it manually.")
        return run_command([conda_path, "env", "create", "-f", env_yaml], f"Creating {env_name} from {env_yaml}")

    return run_command([conda_path, "env", "create", "-f", env_yaml], f"Creating {env_name} from {env_yaml}")

def check_repex_tarean_environment():
    """Optional check for repex_tarean env file.

    Note: REAPER uses the repository-pinned conda env YAMLs under envs/ for workflow execution.
    The upstream repex_tarean repository may (or may not) ship its own environment.yml; it is not
    required for this installer.
    """
    print_section("Checking Repex Tarean Environment")
    
    repex_env_file = "repex_tarean/environment.yml"
    if os.path.exists(repex_env_file):
        print(f"{INFO} Upstream repex_tarean provides: {repex_env_file}")
        return True
    else:
        print(f"{INFO} No repex_tarean/environment.yml found (not required)")
        return False

def install_reportr_environment(conda_path):
    """Install/update the reportr conda environment from envs/reportr.yaml."""
    return _ensure_env_from_yaml(conda_path, "reportr", "envs/reportr.yaml")

def install_repeatexplorer_environment(conda_path):
    """Install/update the repeatexplorer conda environment from envs/repeatexplorer.yaml."""
    return _ensure_env_from_yaml(conda_path, "repeatexplorer", "envs/repeatexplorer.yaml")

def install_reportr_graph_environment(conda_path):
    """Install/update the reportr_graph conda environment from envs/reportr_graph.yaml."""
    return _ensure_env_from_yaml(conda_path, "reportr_graph", "envs/reportr_graph.yaml")

def install_repex_tarean(conda_path):
    """Install repex_tarean tools in the repeatexplorer environment

    Note: Snakemake runs `repex_tarean/seqclust` under the `repeatexplorer` env.
    """
    print_section("Installing Repex Tarean Tools")
    
    # Change to repex_tarean directory
    os.chdir("repex_tarean")
    
    # Compile the source using make within the repeatexplorer environment
    print(f"{INFO} Compiling repex_tarean source code in repeatexplorer environment...")
    success = run_command([conda_path, "run", "-n", "repeatexplorer", "make"], 
                         "Compiling repex_tarean with make", cwd=".")
    
    if not success:
        print(f"{ERROR} Failed to compile repex_tarean")
        os.chdir("..")
        return False
    
    # Check if the compiled binaries exist
    required_binaries = [
        "bin/louvain_community",
        "bin/louvain_convert", 
        "bin/louvain_hierarchy"
    ]
    
    missing_binaries = []
    for binary in required_binaries:
        if not os.path.exists(binary):
            missing_binaries.append(binary)
    
    if missing_binaries:
        print(f"{ERROR} Missing compiled binaries: {missing_binaries}")
        os.chdir("..")
        return False
    
    print(f"{SUCCESS} All required binaries compiled successfully")
    
    # Test the seqclust command within the repeatexplorer environment
    print(f"{INFO} Testing seqclust command...")
    test_success = run_command([conda_path, "run", "-n", "repeatexplorer", "./seqclust", "--help"], 
                              "Testing seqclust help", cwd=".")
    
    if not test_success:
        print(f"{ERROR} seqclust command not working")
        os.chdir("..")
        return False
    
    print(f"{SUCCESS} seqclust command working correctly")
    
    # Optional: run a tiny end-to-end seqclust smoke test if upstream test data exists.
    # The installer should not assume any particular test dataset layout in the cloned repo.
    test_fas = "test_data/LAS_paired_10k.fas"
    if os.path.exists(test_fas):
        print(f"{INFO} Running seqclust smoke test using: {test_fas}")
        os.makedirs("tmp", exist_ok=True)
        test_cmd = [
            conda_path,
            "run",
            "-n",
            "repeatexplorer",
            "./seqclust",
            "-p",
            "-v",
            "tmp/clustering_output",
            test_fas,
        ]
        test_success = run_command(test_cmd, "Running seqclust smoke test", cwd=".")
        if test_success:
            print(f"{SUCCESS} seqclust smoke test passed")
        else:
            print(f"{WARNING} seqclust smoke test failed (skipping; test data/paths may differ upstream)")
    else:
        print(f"{INFO} Skipping seqclust smoke test (missing upstream test file: {test_fas})")
    
    os.chdir("..")
    return True

def create_activation_script(conda_path):
    """Create activation script for easy environment setup"""
    print_section("Creating Activation Script")

    # Write activation helper directly into repo root.
    # This keeps the user-facing entrypoint stable and avoids install_artifacts indirection.
    script_content = f"""#!/usr/bin/env bash
set -euo pipefail

# RepOrtR Environment Activation Script
# Usage:
#   source ./activate_reportr.sh
#   source ./activate_reportr.sh --test

_reportr_echo() {{
  printf '%s\n' "$*"
}}

# This script must be sourced for `conda activate` to affect the current shell.
if [[ "${{BASH_SOURCE[0]}}" == "${{0}}" ]]; then
  _reportr_echo "ERROR: Do not execute this script. Source it instead:"
  _reportr_echo "  source ./activate_reportr.sh"
  exit 1
fi

_reportr_echo "Activating RepOrtR conda environment..."

# Prefer CONDA_EXE if present; otherwise use the installer-discovered path.
CONDA_EXE="${{CONDA_EXE:-{conda_path}}}"

if ! command -v "$CONDA_EXE" >/dev/null 2>&1; then
  if command -v conda >/dev/null 2>&1; then
    CONDA_EXE="conda"
  else
    _reportr_echo "ERROR: conda not found in PATH and CONDA_EXE is not set."
    return 1
  fi
fi

# Initialize conda shell functions if needed (common when sourcing from a non-initialized shell).
if ! declare -F conda >/dev/null 2>&1; then
  _shell_name="${{ZSH_VERSION:+zsh}}"
  if [[ -z "$_shell_name" ]]; then
    _shell_name="bash"
  fi
  eval "$(\"$CONDA_EXE\" \"shell.${{_shell_name}}\" hook)"
fi

conda activate reportr

# Keep seed deterministic for pipeline scripts (Snakemake also sets this from config).
export PYTHONHASHSEED="${{PYTHONHASHSEED:-0}}"

_reportr_echo "✅ RepOrtR environment active"
_reportr_echo "  CONDA_DEFAULT_ENV=${{CONDA_DEFAULT_ENV:-}}"
_reportr_echo "  python=$(command -v python || true)"

if conda env list | awk '{{print $1}}' | grep -qx \"repeatexplorer\"; then
  _reportr_echo "✅ repeatexplorer environment available"
  _reportr_echo "  Run tools in it without switching shells via:"
  _reportr_echo "    conda run -n repeatexplorer <command>"
else
  _reportr_echo "⚠️  repeatexplorer environment not found"
fi

if [[ "${{1:-}}" == "--test" ]]; then
  _reportr_echo "Running REAPER tests (via conda run -n reportr)..."
  if [[ -f "tests/run_tests.py" ]]; then
    "$CONDA_EXE" run -n reportr python tests/run_tests.py
  else
    _reportr_echo "WARNING: tests/run_tests.py not found; skipping."
  fi
fi

_reportr_echo "Ready. Snakemake will select envs via rule-level conda directives."
"""

    script_path = "activate_reportr.sh"
    with open(script_path, "w") as f:
        f.write(script_content)
    
    # Make executable
    os.chmod(script_path, 0o755)

    print(f"{SUCCESS} Activation script created in repo root")
    print(f"{INFO} Use: source activate_reportr.sh")
    print(f"{INFO} Or run tests: source activate_reportr.sh --test")
    
    return True

def verify_installation(conda_path):
    """Verify that all tools are properly installed"""
    print_section("Verifying Installation")
    
    # Test reportr environment
    try:
        result = subprocess.run([conda_path, "run", "-n", "reportr", "python", "-c", 
                               "import Bio, yaml, pandas, numpy, matplotlib; import xlsxwriter; import docx; print('Python packages OK')"],
                              capture_output=True, text=True, timeout=30)
        if result.returncode == 0:
            print(f"{SUCCESS} RepOrtR Python packages verified")
        else:
            print(f"{ERROR} RepOrtR Python packages verification failed")
            return False
    except Exception as e:
        print(f"{ERROR} Could not verify RepOrtR Python packages: {e}")
        return False
    
    # Test repeatexplorer environment
    try:
        result = subprocess.run([conda_path, "run", "-n", "repeatexplorer", "python", "-c", 
                               "print('Repeatexplorer Python OK')"], 
                              capture_output=True, text=True, timeout=30)
        if result.returncode == 0:
            print(f"{SUCCESS} Repeatexplorer environment verified")
        else:
            print(f"{WARNING} Repeatexplorer environment verification failed")
    except Exception as e:
        print(f"{WARNING} Could not verify repeatexplorer environment: {e}")

    # Test reportr_graph environment (R stack for graph reporting)
    try:
        result = subprocess.run(
            [
                conda_path,
                "run",
                "-n",
                "reportr_graph",
                "R",
                "-q",
                "-e",
                'library(igraph); library(visNetwork); cat("reportr_graph R packages OK\\n")',
            ],
            capture_output=True,
            text=True,
            timeout=60,
        )
        if result.returncode == 0:
            print(f"{SUCCESS} reportr_graph environment verified (igraph, visNetwork)")
        else:
            print(f"{WARNING} reportr_graph environment verification failed")
            if result.stderr:
                print(f"{WARNING} R stderr: {result.stderr.strip()}")
    except Exception as e:
        print(f"{WARNING} Could not verify reportr_graph environment: {e}")
    
    # Test bioinformatics tools in reportr environment
    tools_to_test = [
        ("fastqc", ["fastqc", "--version"]),
        ("bbduk", ["bbduk.sh", "version"]),
        ("snakemake", ["snakemake", "--version"])
    ]
    
    for tool_name, cmd in tools_to_test:
        try:
            result = subprocess.run([conda_path, "run", "-n", "reportr"] + cmd, 
                                  capture_output=True, text=True, timeout=30)
            if result.returncode == 0:
                print(f"{SUCCESS} {tool_name} verified in reportr")
            else:
                print(f"{WARNING} {tool_name} not working properly in reportr")
        except Exception as e:
            print(f"{WARNING} Could not verify {tool_name}: {e}")

    # Test seqclust in repeatexplorer env (built in repo under repex_tarean/)
    if os.path.exists("repex_tarean/seqclust"):
        try:
            result = subprocess.run(
                [conda_path, "run", "-n", "repeatexplorer", "repex_tarean/seqclust", "--help"],
                capture_output=True,
                text=True,
                timeout=30,
            )
            if result.returncode == 0:
                print(f"{SUCCESS} seqclust verified (repex_tarean/seqclust)")
            else:
                print(f"{WARNING} seqclust not working properly (repex_tarean/seqclust)")
        except Exception as e:
            print(f"{WARNING} Could not verify seqclust: {e}")
    
    return True

def main():
    """Main installation function"""
    print_header("RepOrtR Installation (Snakemake-managed conda envs)")
    print(f"{INFO} This will ensure three conda envs from repository YAMLs:")
    print(f"   - reportr: workflow + post-processing Python/Snakemake tooling")
    print(f"   - repeatexplorer: RepeatExplorer/TAREAN runtime to build/run repex_tarean/seqclust")
    print(f"   - reportr_graph: R stack for graph reporting (igraph/visNetwork)")
    
    # Check prerequisites
    conda_path = check_conda_available()
    if not conda_path:
        print(f"{ERROR} Conda is required but not available")
        print(f"   Please install conda and run this script again")
        return False
    
    # Clone repex_tarean first if not exists
    if not os.path.exists("repex_tarean"):
        print_section("Cloning Repex Tarean Repository")
        success = run_command(["git", "clone", "https://github.com/kavonrtep/repex_tarean.git"], 
                            "Cloning repex_tarean repository")
        if not success:
            print(f"{ERROR} Failed to clone repex_tarean repository")
            return False
    else:
        print(f"{INFO} repex_tarean directory already exists")
    
    # Optional informational check; does not gate installation.
    check_repex_tarean_environment()
    
    # Install reportr environment
    if not install_reportr_environment(conda_path):
        print(f"{ERROR} Failed to install reportr environment")
        return False
    
    # Install repeatexplorer environment
    if not install_repeatexplorer_environment(conda_path):
        print(f"{ERROR} Failed to install repeatexplorer environment")
        return False

    # Install reportr_graph environment
    if not install_reportr_graph_environment(conda_path):
        print(f"{ERROR} Failed to install reportr_graph environment")
        return False
    
    # Install repex_tarean tools
    if not install_repex_tarean(conda_path):
        print(f"{ERROR} Failed to install repex_tarean tools")
        return False
    
    # Verify installation
    if not verify_installation(conda_path):
        print(f"{WARNING} Some components may not be properly installed")
    
    # Create activation script
    create_activation_script(conda_path)
    
    # Print summary
    print_header("Installation Summary")
    print(f"{SUCCESS} Three environments ensured:")
    print(f"   - reportr (from envs/reportr.yaml)")
    print(f"   - repeatexplorer (from envs/repeatexplorer.yaml)")
    print(f"   - reportr_graph (from envs/reportr_graph.yaml)")
    print(f"{SUCCESS} repex_tarean built: repex_tarean/seqclust")
    print(f"{INFO} Activate with: conda activate reportr")
    print(f"{INFO} Or use: source activate_reportr.sh")
    print(f"{INFO} Run tests: conda run -n reportr python tests/run_tests.py")
    print(f"{INFO} Test seqclust: conda run -n repeatexplorer repex_tarean/seqclust --help")
    print(f"{INFO} Test graph env: conda run -n reportr_graph R -q -e 'library(igraph); library(visNetwork)'")
    
    print(f"\n{INFO} Next steps:")
    print(f"  1. Activate environment: conda activate reportr")
    print(f"  2. Run tests: conda run -n reportr python tests/run_tests.py")
    print(f"  3. Test seqclust: conda run -n repeatexplorer repex_tarean/seqclust --help")
    print(f"  4. Start using RepOrtR pipeline")
    
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 
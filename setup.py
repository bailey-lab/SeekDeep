#!/usr/bin/env python3
import os, subprocess, sys

projectDir = os.path.dirname(os.path.abspath(__file__))
setupToolsDir = os.path.join(projectDir, "cppSetupTools")
if not os.path.exists(os.path.join(setupToolsDir, "setup.py")):
    if not os.path.exists(os.path.join(projectDir, ".git")):
        sys.exit("cppSetupTools is missing; it is a git submodule, so get this project with: git clone --recursive <url>")
    subprocess.run(["git", "submodule", "update", "--init", "cppSetupTools"], cwd=projectDir, check=True)
sys.path.insert(0, setupToolsDir)
from setup import SetupRunner

if __name__ == "__main__":
    SetupRunner.runSetup()

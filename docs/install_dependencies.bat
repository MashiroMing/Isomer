@echo off
REM Install dependencies for MolGenPlus Project

echo ========================================
echo MolGenPlus 依赖安装工具
echo ========================================
echo.

REM Check Python version
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo Error: Python not found. Please install Python 3.8+
    pause
    exit /b 1
)

echo Checking Python installation...
python --version
echo.

echo Please choose installation method:
echo.
echo 1. Conda (Recommended for rdkit) - Best for conda users
echo 2. Pip only - Simpler but may have issues with rdkit
echo 3. Install individual packages manually
echo 0. Cancel
echo.

set /p choice="Please enter your choice (0-3): "

if "%choice%"=="0" (
    echo Installation cancelled.
    pause
    exit /b 0
)

if "%choice%"=="1" (
    echo.
    echo Installing with conda (recommended)...
    echo.

    REM Check if conda is available
    conda --version >nul 2>&1
    if %errorlevel% neq 0 (
        echo Warning: conda not found. Please install Anaconda or Miniconda first.
        echo.
        echo Switching to pip-only installation...
        goto :pip_install
    )

    echo Installing rdkit via conda-forge...
    conda install -c conda-forge rdkit -y
    echo.

    echo Installing Python packages via pip...
    pip install openai pillow scikit-learn numpy pandas

    goto :finish
)

if "%choice%"=="2" (
    echo.
    echo Installing with pip only...
    echo.

    :pip_install
    pip install rdkit openai pillow scikit-learn numpy pandas

    goto :finish
)

if "%choice%"=="3" (
    echo.
    echo Installing individual packages...
    echo.

    set /p packages="Enter packages to install (space-separated): "
    pip install %packages%

    goto :finish
)

echo Invalid choice.
pause
exit /b 1

:finish
echo.
echo ========================================
echo Installation completed!
echo ========================================
echo.
echo Please run 'check_dependencies.py' to verify installation.
echo.
pause

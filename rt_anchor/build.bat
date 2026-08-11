@echo off
REM =====================================================================
REM  build.bat  --  Build the rt-anchor Windows GUI as a single .exe.
REM  Run from the repo root (the folder containing pyproject.toml).
REM  Produces:  dist\rt-anchor.exe
REM
REM  Uses a clean PIP venv on purpose: conda/conda-forge numpy+scipy link
REM  MKL, which roughly doubles the exe size and can throw
REM  "Intel MKL FATAL ERROR" at runtime under onefile extraction.
REM =====================================================================
setlocal enabledelayedexpansion

REM --- Always run from this script's own directory (the repo root) ------
cd /d "%~dp0"

REM --- 1. Locate a Python 3.11 interpreter (needs tkinter) --------------
set "PYEXE="
py -3.11 --version >nul 2>nul && set "PYEXE=py -3.11"
if not defined PYEXE (
    python --version >nul 2>nul && set "PYEXE=python"
)
if not defined PYEXE (
    echo No suitable Python found. Install Python 3.11 (with tkinter^) first.
    exit /b 1
)
echo Using interpreter: %PYEXE%

REM --- 2. Fresh pip-only venv ------------------------------------------
if exist build_venv rmdir /s /q build_venv
%PYEXE% -m venv build_venv || goto :err
call build_venv\Scripts\activate.bat || goto :err

REM --- 3. Tooling + the package with the report extras (matplotlib+plotly)
python -m pip install --upgrade pip wheel setuptools || goto :err
python -m pip install ".[report,app]" || goto :err
python -m pip install "pyinstaller>=6.3" || goto :err

REM --- 4. Guarantee the heavy/unused optional deps are NOT in the env,
REM        so they can never be pulled into the bundle. -----------------
python -m pip uninstall -y rdkit rdkit-pypi kaleido >nul 2>nul

REM --- 5. Clean previous artefacts and build from the spec -------------
if exist build rmdir /s /q build
if exist dist  rmdir /s /q dist
pyinstaller --noconfirm --clean rt_anchor_gui.spec || goto :err

REM --- 6. Drop the .NET app config next to the exe (loadFromRemoteSources) ------
if exist packaging\rt-anchor.exe.config copy /y packaging\rt-anchor.exe.config dist\rt-anchor\rt-anchor.exe.config >nul

echo.
echo ============================================================
echo  BUILD OK  -^>  dist\rt-anchor\rt-anchor.exe  (onedir; zip the folder to ship)
echo  Smoke-test on a machine WITHOUT Python before shipping:
echo    - open a feature table + standards run, Calibrate
echo    - confirm _report.html is several MB and its plots render
echo    - confirm _report.pdf renders text
echo ============================================================
goto :done

:err
echo.
echo  BUILD FAILED (errorlevel %errorlevel%^)
exit /b 1

:done
endlocal

@echo off
REM =====================================================================
REM  build_win.bat -- Build "RT Anchor" desktop as a Windows onedir app.
REM  Lives in build\windows\. Run directly or from the app root; produces:
REM    dist\win\"RT Anchor.exe" (+ _internal\)  and  dist\win\RT Anchor-win64.zip
REM
REM  The release CI (.github\workflows\build-windows.yml) runs this script.
REM
REM  Builds in a pip-only venv at rt_anchor\build_venv, created on first use
REM  with rt-anchor + pywebview + PyInstaller (+ rdkit when it installs).
REM  pip-only on purpose: conda/conda-forge numpy+scipy link MKL, which
REM  roughly doubles the bundle and can throw "Intel MKL FATAL ERROR" frozen.
REM  Delete that folder to rebuild the venv, e.g. after an engine change.
REM =====================================================================
setlocal enabledelayedexpansion

REM --- app root is two levels up (build\windows\) -----------------------
set "ROOT=%~dp0..\.."
set "ENGINE=%~dp0..\..\..\rt_anchor"
set "VENV=%ENGINE%\build_venv"
set "PYEXE=%VENV%\Scripts\python.exe"

REM --- 1. The build venv: create it on first use -----------------------
if not exist "%PYEXE%" (
    set "BASEPY="
    py -3.11 --version >nul 2>nul && set "BASEPY=py -3.11"
    if not defined BASEPY ( python --version >nul 2>nul && set "BASEPY=python" )
    if not defined BASEPY (
        echo No Python found. Install Python 3.11 first.
        exit /b 1
    )
    echo Creating the build venv at !VENV! with !BASEPY!
    !BASEPY! -m venv "%VENV%" || goto :err
    "%PYEXE%" -m pip install --upgrade pip wheel setuptools || goto :err
    pushd "%ENGINE%" || goto :err
    "%PYEXE%" -m pip install ".[report,app]" "pyinstaller>=6.3"
    set "PIPERR=!errorlevel!"
    popd
    if not "!PIPERR!"=="0" goto :err
    "%PYEXE%" -m pip install rdkit || echo WARNING: rdkit did not install - structure hover disabled.
)
echo Using interpreter: %PYEXE%

REM --- 2. Optional rdkit check (molecular-structure hover) --------------
"%PYEXE%" -c "import rdkit" >nul 2>nul
if errorlevel 1 (
    echo WARNING: rdkit not installed in the build env - structure hover disabled.
)

REM --- 3. Clean previous artefacts --------------------------------------
if exist "%ROOT%\dist\win\RT Anchor" rmdir /s /q "%ROOT%\dist\win\RT Anchor"
if exist "%ROOT%\dist\win\RT Anchor-win64.zip" del /q "%ROOT%\dist\win\RT Anchor-win64.zip"
if exist "%~dp0work" rmdir /s /q "%~dp0work"

REM --- 4. Build from the spec (next to this script) ---------------------
"%PYEXE%" -m PyInstaller --noconfirm --clean ^
    --distpath "%ROOT%\dist\win" --workpath "%~dp0work" ^
    "%~dp0RTAnchor.win.spec" || goto :err

REM --- 5. .NET config next to the exe (loadFromRemoteSources fallback) --
copy /y "%~dp0RT Anchor.exe.config" "%ROOT%\dist\win\RT Anchor\RT Anchor.exe.config" >nul

REM --- 6. Shipping zip ---------------------------------------------------
powershell -NoProfile -Command "Compress-Archive -Path '%ROOT%\dist\win\RT Anchor' -DestinationPath '%ROOT%\dist\win\RT Anchor-win64.zip'" || goto :err

echo.
echo ============================================================
echo  BUILD OK
echo    dist\win\"RT Anchor.exe"
echo    dist\win\"RT Anchor-win64.zip"
echo  Smoke-test a windowed ^(no-console^) build by exit code:
echo    "dist\win\RT Anchor\RT Anchor.exe" --selftest
echo  Requirements on target machines: WebView2 Runtime
echo  ^(preinstalled on Win 11 / most Win 10^).
echo ============================================================
goto :done

:err
echo.
echo  BUILD FAILED (errorlevel %errorlevel%^)
exit /b 1

:done
endlocal

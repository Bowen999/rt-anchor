@echo off
REM =====================================================================
REM  build_win.bat -- Build "RT Anchor" desktop as a Windows onedir app.
REM  Lives in build\windows\. Run directly or from the app root; produces:
REM    dist\win\"RT Anchor.exe" (+ _internal\)  and  dist\win\RT Anchor-win64.zip
REM
REM  Requires the rt_anchor build venv (rt_anchor\build_venv) where
REM  rt-anchor + pywebview + PyInstaller + rdkit are pip-installed — see
REM  rt_anchor\build.bat for how that venv is created.
REM =====================================================================
setlocal enabledelayedexpansion

REM --- app root is two levels up (build\windows\) -----------------------
set "ROOT=%~dp0..\.."

REM --- 1. Locate the build venv python ---------------------------------
set "PYEXE=%~dp0..\..\..\rt_anchor\build_venv\Scripts\python.exe"
if not exist "%PYEXE%" (
    echo build_venv not found at %PYEXE%
    echo Create it first: see rt_anchor\build.bat steps 2-3.
    exit /b 1
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

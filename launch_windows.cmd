@echo off
setlocal
pushd "%~dp0"
if errorlevel 1 exit /b 1

set "ACMG_START=import sys;sys.path.insert(0,'.');"
set "ACMG_START=%ACMG_START%from acmg_snv.cli import main;"
set "ACMG_START=%ACMG_START%a=sys.argv[1:];"
set "ACMG_START=%ACMG_START%web=not a or a[0].startswith('-');"
set "ACMG_START=%ACMG_START%a=['web',*a] if web else a;"
set "ACMG_START=%ACMG_START%raise SystemExit(main(a))"
if exist "runtime\python.exe" (
    "runtime\python.exe" -c "%ACMG_START%" %*
    goto finished
)
if exist "python\python.exe" (
    "python\python.exe" -c "%ACMG_START%" %*
    goto finished
)

where py > nul 2> nul
if not errorlevel 1 (
    py -3 -c "%ACMG_START%" %*
    goto finished
)
where python > nul 2> nul
if not errorlevel 1 (
    python -c "%ACMG_START%" %*
    goto finished
)

echo Python 3 is required. Use runtime\python.exe or an installed Python.
pause
popd
exit /b 1

:finished
set "ACMG_EXIT=%ERRORLEVEL%"
if not "%ACMG_EXIT%"=="0" pause
popd
exit /b %ACMG_EXIT%

@echo off
setlocal
chcp 65001 >NUL
set PYTHONUTF8=1
set "JAVA_TOOL_OPTIONS=-Djava.awt.headless=true"
set "PROJECT_DIR=%~dp0\.."
set "PATH=%PROJECT_DIR%\Tool\FastQC;%PATH%"
cd /d "%~dp0"

rem Примеры:
rem   qc_trim.bat                 -> QC без тримминга
rem   qc_trim.bat samplePrefix    -> QC только для образцов с указанным префиксом
rem   qc_trim.bat --trim on       -> QC + тримминг для всех
rem   qc_trim.bat sample --trim on --keep-raw yes

python qc_trim.py %*
pause

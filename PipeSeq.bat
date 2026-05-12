@echo off
setlocal

rem === корень проекта ===
set "PROJECT_DIR=%~dp0"
if "%PROJECT_DIR:~-1%"=="\" set "PROJECT_DIR=%PROJECT_DIR:~0,-1%"

rem === консоль и кодировка ===
chcp 65001 >NUL
set PYTHONUTF8=1
set "JAVA_TOOL_OPTIONS=-Djava.awt.headless=true"

rem === локальный FastQC в PATH (приоритет над системным) ===
set "PATH=%PROJECT_DIR%\Tool\FastQC;%PATH%"

rem === запуск GUI ===
cd /d "%PROJECT_DIR%\Script"
echo Project directory: %PROJECT_DIR%
echo BAT name: %~nx0

python PipeSeq.py

pause

@echo off
REM Получаем путь к текущему BAT-файлу (папка проекта)
set "PROJECT_DIR=%~dp0"
REM Удаляем завершающий обратный слеш, если он есть
if "%PROJECT_DIR:~-1%"=="\" set "PROJECT_DIR=%PROJECT_DIR:~0,-1%"
set "BAT_NAME=%~nx0"

REM Переходим в подпапку Script (где лежат create_shortcut.py и PipeSeq.py)
cd /d "%PROJECT_DIR%\Script"

REM Для отладки можно вывести значения переменных:
echo Project directory: %PROJECT_DIR%
echo BAT name: %BAT_NAME%

REM Запускаем Python-скрипт для создания ярлыка,
REM передавая родительскую папку и имя BAT-файла.
python create_shortcut.py "%PROJECT_DIR%" "%BAT_NAME%"

REM Теперь запускаем основной скрипт PipeSeq.py
python PipeSeq.py

pause

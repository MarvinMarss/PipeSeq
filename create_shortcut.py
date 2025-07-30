import sys
import os
from win32com.client import Dispatch

def main():
    if len(sys.argv) < 3:
        print("Usage: create_shortcut.py <project_directory> <bat_filename>")
        sys.exit(1)
        

    project_dir = sys.argv[1]
    bat_filename = sys.argv[2]
    

    bat_path = os.path.join(project_dir, bat_filename)
    

    icon_path = os.path.join(project_dir, "image.ico")
    

    shortcut_path = os.path.join(project_dir, os.path.splitext(bat_filename)[0] + ".lnk")
    
    shell = Dispatch('WScript.Shell')
    shortcut = shell.CreateShortcut(shortcut_path)
    shortcut.TargetPath = bat_path
    shortcut.WorkingDirectory = project_dir
    if os.path.exists(icon_path):
        shortcut.IconLocation = icon_path
    else:
        print("Файл image.ico не найден, ярлык будет создан без кастомной иконки.")
    shortcut.save()
    print("Ярлык создан:", shortcut_path)

if __name__ == "__main__":
    main()

@echo off
chcp 65001 >nul
echo ========================================
echo   Isomer 主程序打包脚本
echo ========================================
echo.
echo 注意: 打包过程可能需要 5-10 分钟，请耐心等待...
echo.

cd /d "%~dp0"

echo [1/3] 正在清理旧构建文件...
if exist "build" rmdir /s /q "build"
if exist "dist" rmdir /s /q "dist"

echo [2/3] 正在打包 Isomer 主程序...
echo.
python -m PyInstaller Isomer.spec --log-level=WARN
set BUILD_RESULT=%ERRORLEVEL%

echo.
if %BUILD_RESULT%==0 (
    echo [3/3] 打包完成!
    echo.
    echo ========================================
    echo   Isomer 主程序打包成功！
    echo ========================================
    echo.
    echo 可执行文件位于: dist\Isomer\Isomer.exe
    echo.
    echo 是否启动程序进行测试? (Y/N)
    set /p choice=
    if /i "%choice%"=="Y" dist\Isomer\Isomer.exe
) else (
    echo.
    echo ========================================
    echo   打包失败 (错误代码: %BUILD_RESULT%)
    echo ========================================
    echo.
    echo 请检查上面的错误信息
)

echo.
pause

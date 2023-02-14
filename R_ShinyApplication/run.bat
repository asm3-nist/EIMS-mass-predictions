SET ROPTS=--no-save --no-environ --no-init-file --no-restore --no-Rconsole
start http://127.0.0.1:7777
R-Portable\App\R-Portable\bin\Rscript.exe %ROPTS% RunShinyApp.R 1> ShinyApp.log 2>&1

## this portable app building process followed https://www.r-bloggers.com/deploying-desktop-apps-with-r/
message('library paths:\n', paste('... ', .libPaths(), sep='', collapse='\n'))
shiny::runApp('./shiny/wrk/',port=7777)


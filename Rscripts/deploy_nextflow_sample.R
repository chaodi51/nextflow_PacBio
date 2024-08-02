library(optparse)
library(rsconnect)

option_list <- list(
    make_option(c("-r", "--run")),
    make_option(c("-t", "--rmdfile"))
)
opt <- parse_args(OptionParser(option_list=option_list))

appname = opt$run

rsconnect::addConnectServer(url="https://connect.sparkds.io", name="connect", quiet=TRUE)
rsconnect::connectApiUser(account=Sys.getenv("API_USER"), server="connect", apiKey=Sys.getenv("API_KEY"))
rsconnect::deployDoc(doc=opt$rmdfile, appName=appname, account=Sys.getenv("API_USER"), server="connect", launch.browser=FALSE)

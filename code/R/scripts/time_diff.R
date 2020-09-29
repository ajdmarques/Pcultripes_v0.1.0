dat <- read.table(textConnection("start.date start.time end.date end.time
2019-10-11   15:27:26 2019-10-12 06:29:57"), header=TRUE) 

as.numeric(difftime(strptime(paste(dat[,1],dat[,2]),"%Y-%m-%d %H:%M:%S"),
                    strptime(paste(dat[,3],dat[,4]),"%Y-%m-%d %H:%M:%S"),
                    units = 'mins')) 

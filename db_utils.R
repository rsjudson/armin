#--------------------------------------------------------------------------------------
#
# db_utils.R utilities for running MySQL queries
#
# Richard Judson
#
# US EPA
# Questions, comments to: judson.richard@epa.gov
#
# SERVER <<- "mysql-res1.epa.gov"
# USER <<- "rjudson"
#
#--------------------------------------------------------------------------------------
library(RMySQL)
library(DBI)
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
setDBconn <- function(server="mysql-res1.epa.gov",user="rjudson",password) {
  # Sets global variables for the database connection parameters
  #
  # Args:
  #   server: the name of the server (e.g. au.epa.gov)
  #   user: the username
  #   password: the user's password
  #
  # Returns:
  #   no values returned
  #
  printCurrentFunction()
  DB.SERVER <<- server
  DB.USER <<- user
  DB.PASSWORD <<- password
}  
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
getDBconn <- function() {
  # prints global variables for the database connection parameters
  #
  # Args:
  #
  # Returns:
  #   no values returned
  #
  printCurrentFunction()
  if(!exists("DB.SERVER")) cat("DB.SERVER not defined\n")
  else cat("DB.SERVER: ",DB.SERVER,"\n")
  if(!exists("DB.USER")) cat("DB.USER not defined\n")
  else cat("DB.USER: ",DB.USER,"\n")
  if(!exists("DB.PASSWORD")) cat("DB.PASSWORD not defined\n")
  else cat("DB.PASSWORD: ",DB.PASSWORD,"\n")
}  
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
runQuery <- function(query,db,do.halt=F,verbose=F) {
  # Runs a database query and returns a result set
  #
  # Args:
  #   query: a properly formatted SQL query as a string
  #   db: the name of the database
  #   do.halt: if TRUE, halt on errors or warnings
  #   verbose: if TRUE, print diagnostic informaiton
  #
  # Returns:
  #   a result set as a data table
  #
  if(!exists("DB.SERVER")) {
    cat("DB.SERVER not defined\n")
    return(NULL)
  }
  if(!exists("DB.USER")) {
    cat("DB.USER not defined\n")
    return(NULL)
  }
  if(!exists("DB.PASSWORD")) {
    cat("DB.PASSWORD not defined\n")
    return(NULL)
  }
  if(verbose) {
    printCurrentFunction()
    cat("query: ",query,"\n")
    cat("db: ",db,"\n")
  }
  tryCatch({
    con <- dbConnect(drv=RMySQL::MySQL(),user=DB.USER,password=DB.PASSWORD,host=DB.SERVER,dbname=db)
    rs <- dbSendQuery(con, query)
    d1 <- dbFetch(rs, n = -1) 
    if(verbose) {
      print(d1)
      flush.console()
    }
    dbHasCompleted(rs)
    dbClearResult(rs)
    dbDisconnect(con)
    return(d1)
  }, warning = function(w) {
    cat("WARNING:",query,"\n")
    dbDisconnect(con)
    if(do.halt) browser()
    return(NULL)
  }, error = function(e) {
    cat("ERROR:",query,"\n")
    dbDisconnect(con)
    if(do.halt) browser()
    return(NULL)
  })
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
runInsert <- function(query,db,do.halt=F,verbose=F,auto.increment.id=F) {
  # Insert a record into a database
  #
  # Args:
  #   query: a properly formatted SQL query as a string
  #   db: the name of the database
  #   do.halt: if TRUE, halt on errors or warnings
  #   verbose: if TRUE, print diagnostic informaiton
  #   auto.increment: if TRUE, add the auto increment primary key even if not part of the query
  #
  # Returns:
  #   id: if auto.increment=TRUE, return the auto incremented primary key of the record.
  #       Ohterwise, return -1
  #
  if(!exists("DB.SERVER")) {
    cat("DB.SERVER not defined\n")
    return(NULL)
  }
  if(!exists("DB.USER")) {
    cat("DB.USER not defined\n")
    return(NULL)
  }
  if(!exists("DB.PASSWORD")) {
    cat("DB.PASSWORD not defined\n")
    return(NULL)
  }
  if(verbose) {
    printCurrentFunction()
    cat("query: ",query,"\n")
    cat("db: ",db,"\n")
  }
  id <- -1
  tryCatch({
    con <- dbConnect(drv=RMySQL::MySQL(),user=DB.USER,password=DB.PASSWORD,host=DB.SERVER,dbname=db)
    rs <- dbSendQuery(con, query)
    dbHasCompleted(rs)
    dbClearResult(rs)
    if(auto.increment.id) {
      rs2 <- dbSendQuery(con, "select LAST_INSERT_ID()")
      d2 <- dbFetch(rs2, n = -1) 
      id <- d2[1,1]
      dbHasCompleted(rs2)
      dbClearResult(rs2)
    }     
    dbDisconnect(con)
  }, warning = function(w) {
    cat("WARNING:",query," : [",db,"]\n",sep="")
    if(do.halt) browser()
    dbDisconnect(con)
    if(auto.increment.id) return(-1)
  }, error = function(e) {
    cat("ERROR:",query," : [",db,"]\n",sep="")
    print(e)
    if(do.halt) browser()
    dbDisconnect(con)
    if(auto.increment.id) return(-1)
  })
  if(auto.increment.id) return(id)
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
run.insert.table <- function(mat,table,db,do.halt=F,verbose=F) {
  # Inserts multiple rows inot a database table
  #
  # Args:
  #   mat:  data frame containing the data, with the column names corresponding
  #         to the column names of hte database table
  #   table: name of the database table to be inserted in to
  #   db: the name of the database
  #   do.halt: if TRUE, halt on errors or warnings
  #   verbose: if TRUE, print diagnostic informaiton
  #
  # Returns:
  #   no information returned
  #
  if(!exists("DB.SERVER")) {
    cat("DB.SERVER not defined\n")
    return(NULL)
  }
  if(!exists("DB.USER")) {
    cat("DB.USER not defined\n")
    return(NULL)
  }
  if(!exists("DB.PASSWORD")) {
    cat("DB.PASSWORD not defined\n")
    return(NULL)
  }
  if(verbose) {
    printCurrentFunction()
    cat("mat: ",dim(mat),"\n")
    cat("table: ",table,"\n")
    cat("db: ",db,"\n")
  }
  tryCatch({
    con <- dbConnect(drv=RMySQL::MySQL(),user=DB.USER,password=DB.PASSWORD,host=DB.SERVER,dbname=db)
    dbWriteTable(con,name=table,value=mat,field.type=NULL,row.names=F,overwrite=F,append=T)
    dbDisconnect(con)
  }, warning = function(w) {
    cat("WARNING:",table," : [",db,"]\n",sep="")
    dbDisconnect(con)
    if(do.halt) browser()
  }, error = function(e) {
    cat("ERROR:",table," : [",db,"]\n",sep="")
    print(e)
    dbDisconnect(con)
    if(do.halt) browser()
  })
}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
testConnection <- function(db) {
  # Tests the connection to a database
  #
  # Args:
  #   db: the name of the database
  #
  # Returns:
  #   no value returned, but diagnostic informaiton is printed
  #
  printCurrentFunction()
  ret <- runQuery("show databases",db)
  print(ret)
  ret <- runQuery("show tables",db)
  print(ret)
}
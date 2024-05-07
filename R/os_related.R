suppressMessages({
    if (!require(ps)) install.packages("ps")
})


detect_os <- function(){
    os_type <- ps::ps_os_type()
    res <- names(os_type)[os_type]
    return(res)
}


detect_memory <- function(type = c("total", "avail", "percent", "used", "free")){
    type <- match.arg(type)
    res <- ps::ps_system_memory()[[type]]
    return(res)
}


detect_cpu <- function(){
    ps::ps_cpu_count()
}


.os_is_windows <- function(){
    os_name <- Sys.info()[["sysname"]]
    return(os_name)
}

.check_windows_ram <- function(){
    system_info <- system("systeminfo", intern = TRUE)
    phys_mem <- system_info[ grep("Available Physical Memory", system_info) ]
    val <- as.numeric(gsub("\\D*", "", phys_mem, perl = FALSE))
    return(val)
}
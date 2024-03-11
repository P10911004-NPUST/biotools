check_os <- function(){
    os_name <- Sys.info()[["sysname"]]
    return(os_name)
}


check_windows_ram <- function(){
    system_info <- system("systeminfo", intern = TRUE)
    phys_mem <- system_info[ grep("Available Physical Memory", system_info) ]
    val <- as.numeric(gsub("\\D*", "", phys_mem, perl = FALSE))
    return(val)
}

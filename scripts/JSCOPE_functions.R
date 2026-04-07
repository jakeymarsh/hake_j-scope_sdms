## Functions for analysis
## Michael Malick



## add_tranport_id -----------------------------------------
add_transport_id <- function(transport, by = c("lat", "lon1", "lon2"),
                             rm.duplicates = TRUE) {
    ## Add a 'transport_id' column that gives a unique number for each unique
    ## spatial location.
    ##
    ## transport = a data.frame output from `process_transport`
    ## by = vector of column names in `transport` in which to split the data to
    ##      get unique spatial locations
    ## rm.duplicates = logical, if TRUE locations with duplicated transport
    ##                 values are removed, only keeping the first value for the
    ##                 location

    dt <- data.table(transport)
    dt <- dt[ , transport_id:=.GRP, by = by]
    df <- as.data.frame(dt)

    if(rm.duplicates) {
        df.d <- plyr::ddply(df, .(year, ensemble_number, month), transform,
                            dup = duplicated(transport_id))
        df <- df.d[!df.d$dup, ]
        df$dup <- NULL
    }
    row.names(df) <- NULL
    return(df)
}



## process_temperature -------------------------------------
process_temperature <- function(data, initialized, grid, no.se = NULL,
                                fill.bottom.depth = TRUE, ...) {
    ## Process J-SCOPE temperature extractions
    ##
    ## data = data.frame of raw J-SCOPE outputs
    ## initialized = string specifying initialization month
    ## grid = data.frame specifying the J-SCOPE grid
    ## no.se = vector of years with no standard error values, if NULL this is
    ##         ignored
    ## fill.bottom.depth = logical, should locations where target depth < bottom
    ##                     depth be fill with the bottom temperature

    dat <- data

    ## Add read_order column
    dat$read_order <- 1:nrow(dat)

    ## Add initialized column
    dat$initialized <- initialized

    ## lon_rho and lat_rho in april and january initialized forecasts differ at
    ## the ~1e-13 decimal place, presumably due to machine tolerance. This
    ## causes issues later on when joining dataframes. Here, I round the
    ## lon_rho,lat_rho values to avoid joining issues.
    dat <- round_lonlat(dat, digits = 7)
    grd <- round_lonlat(grid, digits = 7)

    ## Join data + grid
    dat.j <- merge(dat, grd, by = c("lat_rho", "lon_rho"),
                   sort = FALSE, all = TRUE)

    ## Apply depth masks
    dat.m <- apply_mask(dat.j)

    ## Only have standard error for 2015 and 2017
    ## Roms output sets the SE to zero in other years, I convert it to NA
    if(!is.null(no.se)) {
        nams <- names(dat.m)[grep("Ste", names(dat.m))]
        for(i in nams) {
            dat.m[[i]] <- ifelse(!dat.m$year %in% no.se, NA, dat.m[[i]])
        }
    }

    ## The ROMS extraction function sets locations where bottom depth < target
    ## depth to NA. Here, I replace those NA values with values from the bottom
    ## most layer.
    if(fill.bottom.depth)
        dat.m <- fill_bottom_depth(dat.m, ...)

    return(dat.m)
}



## round_latlon --------------------------------------------
round_lonlat <- function(data,
                         var.lon = "lon_rho",
                         var.lat = "lat_rho",
                         digits = 12) {
    ## Round the lat,lon columns of input data. Function returns the input
    ## data.frame.
    ## data = a data.frame with lat and lon columns
    ## var.lon = name of lon column in `data`
    ## var.lat = name of lat column in `data`
    ## digits = number of digits to round to

    data[[var.lon]] <- round(data[[var.lon]], digits = digits)
    data[[var.lat]] <- round(data[[var.lat]], digits = digits)
    return(data)
}



## process_transport ---------------------------------------
process_transport <- function(path,
                              pattern = "\\.csv$",
                              ensemble_number = TRUE) {
    ## Process J-SCOPE mass trasport extractions
    ##
    ## path = path to directory with transport data
    ## pattern = regex pattern to match CSV files in `path`
    ## ensemble_number = logical, are the files individual ensemble members
    ##                   These are marked with a "vX" in the file name

    fnames <- list.files(path, pattern = pattern, ignore.case = TRUE)
    out.lst <- vector("list", length(fnames))

    for(i in seq_along(fnames)) {
        i.name <- fnames[i]
        df.i <- read.table(paste0(path, i.name), sep = ",",
                           header = FALSE,
                           stringsAsFactors = FALSE)
        names(df.i) <- c("lat", "lon1", "lon2", "transport_m3_s")
        df.i$read_order  <- 1:nrow(df.i)
        mo.yr            <- filename_year_month(fnames[i],
                                                ensemble_number = ensemble_number)
        df.i$year        <- mo.yr$year
        df.i$month       <- mo.yr$month
        df.i$bimonth     <- mo.yr$bimonth
        df.i$initialized <- "jan"
        if(ensemble_number) df.i$ensemble_number <- mo.yr$ensemble_number
        out.lst[[i]]     <- df.i
    }
    out.df <- plyr::rbind.fill(out.lst)

    ll <- out.df[ , c("lon1", "lon2")]
    out.df$lon_mid <- apply(ll, 1, mean)

    out.df$lat     <- round(out.df$lat, digits = 12)
    out.df$lon1    <- round(out.df$lon1, digits = 12)
    out.df$lon2    <- round(out.df$lon2, digits = 12)
    out.df$lon_mid <- round(out.df$lon_mid, digits = 12)
    return(out.df)
}



## fill_bottom_depth ---------------------------------------
fill_bottom_depth <- function(data,
                              var.fill.anom = c("Temp_Ano_Bot", "Temp_Ano_Bot_v1", "Temp_Ano_Bot_v2", "Temp_Ano_Bot_v3"),
                              var.fill.ste = "Temp_Ste_Bot",
                              vars.100.anom = c("Temp_Ano_100m", "Temp_Ano_100m_v1", "Temp_Ano_100m_v2", "Temp_Ano_100m_v3"),
                              vars.100.ste  = c("Temp_Ste_100m", "Temp_Ste_100to400m", "Temp_Ste_100to400m_Avg"),
                              vars.250.anom = c("Temp_Ano_250m", "Temp_Ano_250m_v1", "Temp_Ano_250m_v2", "Temp_Ano_250m_v3"),
                              vars.250.ste  = c("Temp_Ste_250m")) {
    ## Fill locations where bottom depth < target variable depth with data from
    ## the `var.fill` columns. The function returns the input data frame with
    ## additional variable columns appended with "_fill".
    ##
    ## data = a data.frame
    ## var.fill.anom = column name of data used to fill shallow area anomalies
    ## var.fill.ste = column name of data used to fill shallow areas standard
    ##                errors
    ## vars.100.anom = 100m anomaly column names to fill
    ## vars.100.ste = 100m standard error column names to fill
    ## vars.250.anom = 100m anomaly column names to fill
    ## vars.250.ste = 100m standard error column names to fill

    dat <- data

    if(!is.null(vars.100.anom)) {
        for(i in seq_along(vars.100.anom)) {
            v_fill <- var.fill.anom[i]
            vv <- vars.100.anom[i]
            vv.n <- paste0(vv, "_fill")
            dat[[vv.n]] <- ifelse(dat$exclude_100m, dat[[v_fill]], dat[[vv]])
        }
    }
    if(!is.null(vars.100.ste)) {
        for(i in seq_along(vars.100.ste)) {
            vv <- vars.100.ste[i]
            vv.n <- paste0(vv, "_fill")
            dat[[vv.n]] <- ifelse(dat$exclude_100m, dat[[var.fill.ste]], dat[[vv]])
        }
    }

    if(!is.null(vars.250.anom)) {
        for(i in seq_along(vars.250.anom)) {
            v_fill <- var.fill.anom[i]
            vv <- vars.250.anom[i]
            vv.n <- paste0(vv, "_fill")
            dat[[vv.n]] <- ifelse(dat$exclude_250m, dat[[v_fill]], dat[[vv]])
        }
    }
    if(!is.null(vars.250.ste)) {
        for(i in seq_along(vars.250.ste)) {
            vv <- vars.250.ste[i]
            vv.n <- paste0(vv, "_fill")
            dat[[vv.n]] <- ifelse(dat$exclude_250m, dat[[var.fill.ste]], dat[[vv]])
        }
    }

    return(dat)
}



## apply_mask ----------------------------------------------
apply_mask <- function(data) {
    ## Apply depth mask to vars in data
    ##
    ## This function sets values of a variable to NA if the
    ## bottom depth < target variable depth

    dat <- data

    vars.100 <- names(dat)[grep("^Temp.*100", names(dat))]
    vars.250 <- names(dat)[grep("^Temp.*250", names(dat))]

    for(i in seq_along(vars.100)) {
        vv <- vars.100[i]
        dat[[vv]] <- ifelse(dat$exclude_100m, NA, dat[[vv]])
    }

    for(i in seq_along(vars.250)) {
        vv <- vars.250[i]
        dat[[vv]] <- ifelse(dat$exclude_250m, NA, dat[[vv]])
    }

    return(dat)
}



## write_jscope_vars ---------------------------------------
write_jscope_vars <- function(data, out.dir) {
    ## Write processed JSCOPE data to csv files
    ##
    ## One file is created for each unique entry in `filename` column in data
    ##
    ## data = a data.frame
    ## out.dir = directory to save csv files

    if(!dir.exists(out.dir))
        dir.create(out.dir, recursive = TRUE)

    sp <- split(data, data$filename)

    for(i in 1:length(sp)) {
        sp.i <- sp[[i]]
        f.name <- gsub("nc$", "csv", unique(sp.i$filename))
        f.path <- paste0(out.dir, f.name)
        write.table(sp.i, file = f.path, sep = ",", row.names = FALSE)
    }
}



## read_jscope_vars  ---------------------------------------
read_jscope_vars <- function(path, pattern = "\\.nc$", return.df = FALSE,
                             year_month = TRUE) {
    ## Read JSCOPE re-forecast variables
    ##
    ## This function reads all ncdf files in the dir, converts them to
    ## data.frames, and stores them in a list (each list component cooresponds
    ## to a single file) and returns the list or rbinded data.frame.
    ##
    ## path = path to root directory with netcdf files
    ## pattern = regex for filename to match
    ## return.df = logical, if TRUE the function returns a data.frame, otherwise
    ##             it returns a list
    ## year_month = logical, should the year and month be read from the file
    ##              name

    if(!dir.exists(path))
        stop("`path` doesn't exist", call. = FALSE)

    fls           <- list.files(path, pattern = pattern)
    out.ls        <- vector("list", length(fls))
    fls.noext     <- gsub("\\.nc$", "", fls)
    names(out.ls) <- fls.noext

    for(i in seq_along(fls)) {
        fls.i    <- fls[i]
        lst      <- read_ncdf_list(paste0(path, fls.i))
        lst.vec  <- lapply(lst, as.vector)
        df       <- as.data.frame(lst.vec, stringsAsFactors = FALSE)
        df$read_order <- 1:nrow(df)
        if(year_month) {
            add.columns   <- filename_year_month(lst$filename)
            df            <- cbind(df, add.columns)
        }
        df$filename   <- basename(df$filename)
        df[is.na(df)] <- NA
        out.ls[[i]]   <- df
    }
    cat("Read", length(fls), "files", "\n")

    if(return.df) out.ls <- plyr::rbind.fill(out.ls)
    return(out.ls)
}



## filename_year_month -------------------------------------
filename_year_month <- function(filename, ensemble_number = FALSE) {
    ## Extract the year, month, and bi-monthly indicator for a file name
    ## The function returns a data.frame with columns:
    ##      $year
    ##      $month
    ##      $bimonth
    ##
    ## filename = a filename to parse

    f.split <- filename_split(filename)

    ## extract year
    year <- f.split[grep("[0-9]{4}", f.split)]

    if(length(year) == 0)
        warning("Year not found in file name", call. = FALSE)

    ## extract month
    month.vec <- c("Jun", "Jul", "Aug", "Sep")
    month     <- unique(grep(paste(month.vec, collapse = "|"),
                             f.split, value = TRUE))

    if(length(month) == 0)
        warning("Month not found in file name", call. = FALSE)

    ## bi-monthly
    bi         <- gsub(paste(month.vec, collapse = "|"), "", month)
    month.half <- ifelse(bi == "", NA, bi)
    month      <- ifelse(!is.na(month.half),
                         gsub(paste0(month.half, "$"), "", month),
                         month)

    out.df <- data.frame(year = as.numeric(year),
                         month = month,
                         bimonth = month.half,
                         stringsAsFactors = FALSE)

    ## ensemble number
    if(ensemble_number) {
        en <- grep("^v[0-9]$", f.split, value = TRUE)
        out.df$ensemble_number <- en
    }

    return(out.df)
}



## filename_split ------------------------------------------
filename_split <- function(filename, split = "_") {
    ## Split a file name
    ##
    ## This function splits a filename according to the pattern specified in
    ## `split` and returns the splitted file name as a vector.
    ##
    ## filename = filename to split, can be a path to a file, in which case the
    ##            path component is removed
    ## split = pattern to split the filename by
    fn.base   <- basename(filename)
    fn.noext  <- gsub("\\..*$", "", fn.base)
    fn.vector <- strsplit(fn.noext, "_")
    return(unlist(fn.vector))
}

if(FALSE) {

    fn <- "test_test.nc"
    filename_split(fn)

    fn <- "test_test-test.nc"
    filename_split(fn)

}



## jscope_combine_fore_hind --------------------------------
jscope_combine_fore_hind <- function(forecast, hindcast) {
    ## Combine forecast and hindcast temperature data.frames
    ##
    ## forecast = data.frame output from `jscope.dataframe.temp`
    ## hindcast = data.frame output from `jscope.dataframe.temp`

    out <- forecast
    out$temp_hind <- hindcast$temp
    out$temp <- NULL
    out$temp_fore <- forecast$temp
    out$temp_diff <- forecast$temp - hindcast$temp
    out$filename_fore <- out$filename
    out$filename_hind <- hindcast$filename
    out$filename <- NULL
    return(out)
}

if(FALSE) {

    p.hind <- "./data/jscope-outputs-2017-09-25/forMike/2015_Hindcast_Jun_forMike.nc"
    p.fore <- "./data/jscope-outputs-2017-09-25/forMike/2015_Reforecast_Jun_forMike.nc"

    vars <- c("temp", "lon_rho", "lat_rho", "s_rho", "depth")
    hind <- read_ncdf_list(p.hind, vars = vars)
    fore <- read_ncdf_list(p.fore, vars = vars)


    ## Only `temp` and `filename` should differ
    all.equal(hind, fore)

    ## Convert to data.frame
    hind.df <- jscope_dataframe_3d(hind, 2015, 6)
    fore.df <- jscope_dataframe_3d(fore, 2015, 6)

    xx <- jscope_combine_fore_hind(fore.jun.2015.df, hind.jun.2015.df)
    head(xx)
    tail(xx)
    sapply(xx, class)

}



## jscope_dataframe_3d -------------------------------------
jscope_dataframe_3d <- function(data, year, month,
                                temp.var = "temp",
                                depth.var = "depth",
                                lon.var = "lon_rho",
                                lat.var = "lat_rho") {
    ## Convert 3D JSCOPE data read in using `read_ncdf_list` to a data.frame
    ##
    ## data = list output from `read_ncdf_list()`
    ## year = year of data
    ## month = month of data
    ## temp.var = name of temperature variable
    ## depth.var = name of depth variable
    ## lon.var = name of longitude variable
    ## lat.var = name of latitude variable

    temp  <- data[[temp.var]]
    depth <- data[[depth.var]]
    layers <- dim(temp)[3]

    ## Name dimensions
    temp.names <- list(data[[lon.var]][ , 1],
                       data[[lat.var]][1, ],
                       1:layers)
    depth.names <- list(data[[lon.var]][ , 1],
                        data[[lat.var]][1, ],
                        1:layers)

    dim.nm <- c(lon.var, lat.var, "layer")
    names(temp.names) <- dim.nm
    names(depth.names) <- dim.nm

    dimnames(temp) <- temp.names
    dimnames(depth) <- depth.names


    ## Convert temp to data.frame
    xt <- as.data.frame.table(temp, responseName = temp.var,
                              stringsAsFactors = FALSE)
    xt[[lon.var]] <- as.numeric(xt[[lon.var]])
    xt[[lat.var]] <- as.numeric(xt[[lat.var]])
    xt$layer <- as.numeric(xt$layer)


    ## Convert depth to data.frame
    xd <- as.data.frame.table(depth, responseName = depth.var,
                              stringsAsFactors = FALSE)
    xd[[lon.var]] <- as.numeric(xd[[lon.var]])
    xd[[lat.var]] <- as.numeric(xd[[lat.var]])
    xd$layer <- as.numeric(xd$layer)

    ## Join temp and depth
    xx <- plyr::join(xt, xd, by = dim.nm)

    ## Add unique cell id
    ## xx <- plyr::ddply(xx, .(layer), transform,
    ##                   cell.id = 1:length(lon))

    ## Convert NaN to NA
    xx[[depth.var]][is.nan(xx[[depth.var]])] <- NA

    out <- data.frame(year = year,
                      month = month,
                      layer = xx$layer,
                      lon = xx[[lon.var]],
                      lat = xx[[lat.var]],
                      ## cell_id = xx$cell.id,
                      depth = xx[[depth.var]],
                      temp = xx[[temp.var]],
                      filename = data$filename,
                      stringsAsFactors = FALSE)
    names(out) <- c("year", "month", "layer", lon.var, lat.var, depth.var,
                    temp.var, "filename")
    return(out)
}

if(FALSE) {

    pth <- "./data/jscope-outputs-2017-09-25/forMike/2015_Hindcast_Jun_forMike.nc"
    vars <- c("temp", "lon_rho", "lat_rho", "s_rho", "depth")
    xa <- read_ncdf_list(pth, vars = vars)
    xx <- jscope_dataframe_3d(xa, 2015, 6)
    head(xx)
    tail(xx)
    sapply(xx, class)

}



## read_ncdf_list ------------------------------------------
read_ncdf_list <- function(ncdf, filename = "filename", vars = NULL, ...) {
    ## Read variables in a ncdf file and store in a list
    ##
    ## This function takes as input a path to a netcdf file and
    ## reads variables in the ncdf file and writes them to an R list, which
    ## is returned. The names of the output list are the same as the variable
    ## names in the ncdf file. If `vars` is NULL, all variables in the netcdf
    ## file are read.
    ##
    ## ncdf = path to a netcdf file
    ## filename = if non-null, add the file name to output list where the name
    ##            of the list component is specified by `filename`
    ## vars = char vector specifying the variable names in `ncdf` to read
    ## ... = passed to ncvar_get()

    ## Open connection to netcdf file
    con.ncdf <- ncdf4::nc_open(ncdf)

    if(is.null(vars)) {
        vars <- attributes(con.ncdf$var)$names
    }

    lst <- vector("list", length(vars))
    names(lst) <- vars

    for(i in seq_along(vars)) {
        lst[[i]] <- ncdf4::ncvar_get(con.ncdf, varid = vars[i], ...)
    }

    if(!is.null(filename)) {
        lst[[filename]] <- basename(con.ncdf$filename)
    }

    ## Close connection
    ncdf4::nc_close(con.ncdf)

    return(lst)
}

## Functions for analysis


## raster_grid_var -----------------------------------------
raster_grid_var <- function(data, grid, var, var_new, xvar, yvar, fun, by = NULL, ...) {
    ## Grid data using specified grid.
    ##
    ## This functions takes a predefined raster grid and a data frame of data
    ## and calculates per-grid-cell statistics. The function returns a
    ## data frame of the gridded data set.
    ##
    ## data = a data.frame
    ## grid = a RasterLayer object (e.g., output from grid_create)
    ## var = column name in 'data' of variable to grid
    ## var_new = new variable name
    ## xvar = column name of 'data' giving x-coordinate
    ## yvar = column name of 'data' giving y-coordinate
    ## fun = function to summarize data within a cell, see ?raster
    ## by = character vector of variables in 'data' to split the data by.
    ##      Separate grids are created for each subset of data.
    ## ... = passed to 'rasterize'

    if(is.null(by)) {
        out <- .raster_grid_var(data = data,
                                grid = grid,
                                var = var,
                                var_new = var_new,
                                xvar = xvar,
                                yvar = yvar,
                                fun = fun, ...)
    } else {
        out <- plyr::ddply(data, by, function(x) {
            out <- .raster_grid_var(data = x,
                                    grid = grid,
                                    var = var,
                                    var_new = var_new,
                                    xvar = xvar,
                                    yvar = yvar,
                                    fun = fun, ...)
        })
    }

    return(out)
}

.raster_grid_var <- function(data, grid, var, var_new, xvar, yvar, fun, ...) {

    mat <- matrix(c(data[[xvar]], data[[yvar]]), byrow = FALSE, ncol = 2)
    ras <- raster::rasterize(x = mat, y = grid, field = data[[var]],
                     fun = fun, ...)

    ## Create a unique cell identifier
    grid_df <- as.data.frame(grid, xy = TRUE)
    grid_df[["cell_id"]] <- 1:nrow(grid_df)
    cellid <- raster::rasterize(x = grid_df[ , c("x", "y")],
                        y = grid, field = grid_df[["cell_id"]],
                        fun = function(x, ...) unique(x))

    ## Create raster stack
    lst <- list(cellid, ras)
    names(lst) <- c("cell_id", var_new)
    grd_stack <- raster::stack(lst)

    ## Convert to data frame
    grd_df <- as.data.frame(grd_stack, xy = TRUE)
    names(grd_df)[names(grd_df) == "x"] <- xvar
    names(grd_df)[names(grd_df) == "y"] <- yvar

    return(grd_df)
}


library(raster)
## raster_grid_create --------------------------------------
raster_grid_create <- function(data, res, crs,
                               xvar = "easting_m",
                               yvar = "northing_m",
                               xmin = NULL,
                               xmax = NULL,
                               ymin = NULL,
                               ymax = NULL) {
    ## Create a raster grid. The default is to use the spatial extent of the
    ## input data for the grid extent, but setting xmin,xmax,ymin,ymax overrides
    ## this.
    ##
    ## Returns a 'RasterLayer' object
    ##
    ## data = a data.frame
    ## res = grid resolution in same units as input data x,y
    ## crs = character or object of class CRS, see ?raster
    ## xvar = column name of 'data' giving x-coordinate
    ## yvar = column name of 'data' giving y-coordinate
    ## xmin = min x-coord extent, if NULL min xvar value is used
    ## xmax = max x-coord extent, if NULL max xvar value is used
    ## ymin = min y-coord extent, if NULL min yvar value is used
    ## ymax = max y-coord extent, if NULL max yvar value is used

    xmin <- ifelse(is.null(xmin), min(data[[xvar]]), xmin)
    xmax <- ifelse(is.null(xmax), max(data[[xvar]]), xmax)
    ymin <- ifelse(is.null(ymin), min(data[[yvar]]), ymin)
    ymax <- ifelse(is.null(ymax), max(data[[yvar]]), ymax)

    grd <- raster::raster(xmn = xmin, xmx = xmax,
                          ymn = ymin, ymx = ymax,
                          res = res, crs = crs)

    return(grd)
}



## get_edf -------------------------------------------------
gam_edf <- function(fits, pattern, names = NULL) {
    ## Get EDF and significance for a term in a list of GAM models.

    lst <- lapply(fits, .gam_edf, pattern = pattern) ## list of fits
    out <- as.data.frame(do.call("rbind", lst))
    names(out) <- "edf"

    if(is.null(names(fits))) {
        model <- paste0("M", 0:(length(fits) - 1))
    } else {
        model <- names(fits) ## list of fits
    }
    out <- cbind(model, out)

    if(!is.null(names) & length(names) == nrow(out)) out[ , "model"] <- names

    row.names(out) <- NULL
    return(out)
}

if(FALSE) {

    gam_edf(gam_po, "temp_100")

}

.gam_edf <- function(fit, pattern) {
    ## Get EDF and significance code for a GAM model term. Function returns a
    ## vector of the EDF with a "*" appended if p < 0.05.
    ##
    ## fit = a GAM model object
    ## pattern = pattern to select model term
    if(!is(fit, "gam")) {
        stop("Input is not class 'gam'")
    }

    x <- summary(fit)$s.table
    s <- grepl(pattern, row.names(x))

    edf  <- x[s, "edf"]
    pval <- x[s, "p-value"]
    if(length(edf) > 0) {
        star <- ifelse(pval < 0.05, "*", "")
        out <- paste0(formatC(edf, digits = 2, format = "f"), star)
    } else {
        out <- "-"
    }
    return(out)
}

if(FALSE) {

    .gam_edf(gam_po$po3, "v_100_300_kr")
    .gam_edf(gam_po$po0, "v_100_300_kr")

}


## gam_summarize -------------------------------------------
gam_summarize <- function(fits, order.aic = TRUE, auc = FALSE, names = NULL) {
    ## Summarize GAM model fits
    ##
    ## This function summarizes a fitted GAM model by pulling out relevant
    ## output stored in the GAM object. The function outputs a data.frame with a
    ## single row for each model input. If more than one model is supplied, the
    ## delta-AIC across models is also calculated.
    ## See: ?gamObject
    ## See: getAnywhere(AIC.default)
    ##
    ## fit = a list of fitted GAM model (object of class 'gam')
    ## order.aic = logical, should results be sorted by lowest AIC
    ## auc = logical, should an AUC be calculated
    ## names = optional vector of names to use for models

    lst <- lapply(fits, .gam_summarize, auc = auc) ## list of fits
    out       <- plyr::rbind.fill(lst)
    min.ind   <- which.min(out$aic)
    out$d.aic <- out$aic - out$aic[which.min(out$aic)]
    if(is.null(names(fits))) {
        model <- paste0("M", 0:(length(fits) - 1))
    } else {
        model <- names(fits) ## list of fits
    }
    out <- cbind(model, out)

    if(!is.null(names) & length(names) == nrow(out)) out[ , "model"] <- names

    if(order.aic & length(out$d.aic) > 1) {
        out <- out[order(out$d.aic), ]
    }

    row.names(out) <- NULL
    return(out)
}


.gam_summarize <- function(x, auc = FALSE) {
    ## Summarize GAM model fit
    ##
    ## This function summarizes a fitted GAM model by pulling out relevant
    ## output stored in the GAM object. The function outputs a data.frame with a
    ## single row. See: ?gamObject
    ##
    ## x = a fitted GAM model (object of class 'gam')
    ## auc = logical, should an AUC be calculated

    if(!is(x, "gam")) {
        stop("Input is not class 'gam'")
    }

    aic           <- x$aic
    deviance      <- x$deviance
    null.deviance <- x$null.deviance
    converged     <- x$converged
    boundary      <- x$boundary
    r.sq          <- summary(x)$r.sq
    gcv           <- x$gcv.ubre
    method        <- x$method
    edf           <- sum(summary(x)$edf)
    n             <- summary(x)$n
    np            <- summary(x)$np
    dev.expl      <- summary(x)$dev.expl
    family        <- x$family$family
    link          <- x$family$link
    mse           <- mean(residuals(x, type = "response")^2)
    rmse          <- sqrt(mse)

    edf <- ifelse(length(edf) == 0, NA, edf)

    df <- data.frame(dev.expl      = dev.expl,
                     r.sq          = r.sq,
                     aic           = aic,
                     deviance      = deviance,
                     gcv           = gcv,
                     method        = method,
                     null.deviance = null.deviance,
                     converged     = converged,
                     boundary      = boundary,
                     edf           = edf,
                     n             = n,
                     np            = np,
                     family        = family,
                     link          = link,
                     mse           = mse,
                     rmse          = rmse,
                     stringsAsFactors = FALSE)

    if(auc) {
        df$auc <- pROC::roc(x$y, x$fitted.values, plot = FALSE, quiet = TRUE)$auc[1]
    }

    row.names(df) <- NULL
    return(df)
}


## gam_vc --------------------------------------------------
gam_vc <- function(model, newdata, var, term) {
    ## Get 2D variable coefficient slope for prediction grid
    ##
    ## model = fitted gam model
    ## newdata = newdata used for prediction
    ## var = string, column name in newdata for the vc variable
    ## term = string, gam model term for vc term

    pred <- predict(model, newdata = newdata, type = "terms",
                    terms = term, se = TRUE)

    pdat <- newdata
    pdat$slope <- as.vector(pred$fit[ , 1] / newdata[ , var])
    pdat$ci    <- as.vector((pred$se[ , 1] / newdata[ , var]) * qnorm(0.975))
    pdat$lower <- pdat$slope - pdat$ci
    pdat$upper <- pdat$slope + pdat$ci
    pdat$sig   <- (pdat$lower < 0 & pdat$upper < 0) | (pdat$lower > 0 & pdat$upper > 0)

    xmin <- min(newdata[["easting_m"]], na.rm = TRUE)
    xmax <- max(newdata[["easting_m"]], na.rm = TRUE)
    ymin <- min(newdata[["northing_m"]], na.rm = TRUE)
    ymax <- max(newdata[["northing_m"]], na.rm = TRUE)
    xmin <- xmin - (xmin * 0.01)
    xmax <- xmax + (xmin * 0.01)
    ymin <- ymin - (ymin * 0.01)
    ymax <- ymax + (ymin * 0.01)
    limits = data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)

    list(data = pdat, limits = limits)

}



## plot_gam_2d_vc ------------------------------------------
plot_gam_2d_vc <- function(x,
                           title = "",
                           xlab = "Eastings",
                           ylab = "Northings",
                           fill_lab = "",
                           point_size = 0.7,
                           lonlat = FALSE,
                           plot = TRUE,
                           se = TRUE) {
    ## Plot variable coefficient GAM slopes
    ##
    ## x = list of data, should be in format output from `plot.gam`
    ## type = "bubble" or "point" or "raster"
    ## se = logical, should significant locations be indicated (only implemented
    ##      for type = "raster"

    col.mid   <- "grey95"
    col.low   <- RColorBrewer::brewer.pal(8, "Set1")[2]
    col.high  <- RColorBrewer::brewer.pal(8, "Set1")[1]

    gr       <- expand.grid(x = x$x, y = x$y)
    pdat     <- data.frame(gr, slope = as.vector(x$fit))
    pdat     <- na.omit(pdat)
    pdat$col <- ifelse(pdat$slope > 0, "pos", "neg")

    xmin <- min(pdat[["x"]], na.rm = TRUE)
    xmax <- max(pdat[["x"]], na.rm = TRUE)
    ymin <- min(pdat[["y"]], na.rm = TRUE)
    ymax <- max(pdat[["y"]], na.rm = TRUE)
    xmin <- xmin - (xmin * 0.01)
    xmax <- xmax + (xmin * 0.01)
    ymin <- ymin - (ymin * 0.01)
    ymax <- ymax + (ymin * 0.01)

    if(lonlat) {
        map <- namerica_state_wgs_df
        bath <- bath_200_wgs_df
    } else {
        map <- namerica_state_utm_df
        bath <- bath_200_utm_df
    }

    if(se) {
        ci <- na.omit(x$se * qnorm(0.975))
        pdat$lower <- pdat$slope - ci
        pdat$upper <- pdat$slope + ci
        pdat$sig <- (pdat$lower < 0 & pdat$upper < 0) |
                    (pdat$lower > 0 & pdat$upper > 0)
    }

    g <- ggplot(pdat) +
        geom_tile(aes(x, y, fill = slope))
    if(se) {
        g <- g + geom_point(data = pdat[pdat$sig, ], aes(x = x, y = y),
                            size = point_size, shape = 16)
    }
    g <- g + geom_sf(data = map, fill = "grey90", color = "grey30", size = 0.2) +
        scale_x_continuous(limits = c(xmin, xmax), expand = expansion(mult = 0.05)) +
        scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = 0.05)) +
        ggtitle(title) +
        scale_fill_gradient2(low = col.low,
                             mid = col.mid,
                             high = col.high, midpoint = 0.0,
                             space = "Lab",
                             na.value = "green") +
        labs(fill = fill_lab, x = xlab, y = ylab) +
        theme_simple() +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank())

    if(plot) print(g)
    invisible(pdat)
}



## plot_gam_1d ---------------------------------------------
plot_gam_1d <- function(x, se.mult = 1,
                        xlab = NULL,
                        ylab = NULL,
                        title = NULL,
                        rug = TRUE,
                        size = 1,
                        plot = TRUE) {
    ## Plot a GAM 1-d smooth
    ##
    ## x = list of data, should be in format output from `plot.gam`
    ## se.mult = value to multiply standard error by
    ## rug = logical, should a rug be added to plot

    if("y" %in% names(x)) {
        stop("x is not a 1-d smooth")
    }

    fit    <- as.vector(x$fit)
    se.fit <- x$se
    lower  <- fit - (se.mult * se.fit)
    upper  <- fit + (se.mult * se.fit)

    if(is.null(xlab)) xlab <- x$xlab
    if(is.null(ylab)) ylab <- x$ylab

    dat <- data.frame(fit = as.vector(x$fit),
                      se.fit = x$se,
                      x = x$x,
                      lower = fit - (se.mult * se.fit),
                      upper = fit + (se.mult * se.fit))

    g <- ggplot(dat) +
        geom_hline(yintercept = 0, color = "grey50", linetype = 2, size = 0.3) +
        geom_ribbon(aes(x = x, ymin = lower, ymax = upper), fill = "grey25", alpha = 0.2) +
        geom_line(aes(x = x, y = fit), size = size) +
        labs(x = xlab,
             y = ylab,
             title = title) +
        theme_simple()

    if(rug) {
        g <- g + geom_rug(data = data.frame(x = x$raw), aes(x = x),
                          length = unit(0.02, "npc"),
                          color = "grey40")
    }

    if(plot) print(g)

    invisible(dat)
}



add_shelf_distance <- function(data,
                               shelf,
                               proj.data      = proj_wgs,
                               proj.shelf     = proj_wgs,
                               proj.dist      = proj_utm,
                               var.shelf.dist = "shelf_distance_m",
                               lon.var        = "lon",
                               lat.var        = "lat") {
  ## Calculate distance b/w points and a bathymetric polygon
  ##
  ## This function takes as input a data.frame of point data (i.e., `data`)
  ## and a spatial polygon object (i.e., `shelf`) and calculates the distance
  ## b/w each point and the polygon boundary. Additionally, whether the point
  ## is located inside the polygon ("off") or outside ("on") is recorded. For
  ## bathymetric polygons, "off" indicates the point is off the shelf and "on"
  ## indicates the point is on the shelf.
  ##
  ## The function returns the input `data` data.frame with one additional
  ## column:
  ##   1. var.shelf.dist = distance to shelf in units of `proj.dist`
  ##
  ##      + values = off shelf (seaward of shelf)
  ##      - values = on shelf (shorward of shelf)
  ##
  ## data = data.frame to add distance to shelf to
  ## shelf = sp object of bathymetry
  ## proj.data = projection of input `data`
  ## proj.shelf = projection of input `shelf` data, expected to be `sp` class
  ## proj.dist = projection to use for calculating distance,
  ##             distances can't be calculted in lon,lat units
  ## var.shelf.dist = name of new column in `data` specifying distance to shelf
  ## lon.var = char, column name in `data` specifying longitude
  ## lat.var = char, column name in `data` specifying latitude
  
  if(!lon.var %in% names(data))
    stop("`lon.var` not found")
  if(!lat.var %in% names(data))
    stop("`lat.var` not found")
  
  dat <- data
  
  ## Create internal row indicator
  dat$ind_shelf_internal <- 1:nrow(dat)
  
  ## Project input data
  dat.ll <- subset(dat, select = c(lon.var, lat.var, "ind_shelf_internal"))
  dat.ll <- sf::st_as_sf(dat.ll, coords = c(lon.var, lat.var))
  st_crs(dat.ll) <- proj.data
  dat.proj <- sf::st_transform(dat.ll, proj.dist)
  
  ## Project shelf data
  st_crs(shelf) <- proj.shelf
  shelf.proj <- sf::st_transform(shelf, proj.dist)
  shelf.line <- st_cast(shelf.proj, "MULTILINESTRING")
  
  ## Calc distance to bathy
  dist.all <- sf::st_distance(dat.proj, shelf.line, by_element = FALSE)
  dist.all <- matrix(dist.all, nrow = nrow(dist.all))
  dist <- apply(dist.all, 1, min)
  dat[[var.shelf.dist]] <- as.vector(dist)
  
  ## Inside bathy polygon == off shelf
  shelf.off <- dat.proj[shelf.proj, ]
  dat[[var.shelf.dist]] <- ifelse(dat$ind_shelf_internal %in%
                                    shelf.off$ind_shelf_internal,
                                  dat[[var.shelf.dist]],
                                  dat[[var.shelf.dist]] * -1)
  
  ## Clean-up
  dat$ind_shelf_internal <- NULL
  
  return(dat)
}





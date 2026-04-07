## Functions for hake analysis
## Michael Malick


## surf_plot -----------------------------------------------
surf_plot <- function(surf,
                      main = "",
                      var.x = "x",
                      var.y = "y",
                      var.z = "z",
                      var.year = "year",
                      guide.title = "",
                      guide.lims = NULL,
                      map = FALSE,
                      map.data = namerica_utm_df) {
    ## Plot fitted surface
    ##
    ## surf = output from surf.temp()
    ## main = title for plot
    ## map = logical, overlay a map
    ## map.data = data.frame of map data

    df <- data.frame(x = surf[[var.x]],
                     y = surf[[var.y]],
                     z = surf[[var.z]],
                     year = surf[[var.year]])

    g <- ggplot(df)

    if(map) {
        g <- g + geom_map(data = map.data,
                          map = map.data,
                          aes(map_id = id),
                          fill = "white", color = "grey10", size = 0.2)
    }

    g <- g + geom_raster(aes(x = x, y = y, fill = z)) +
        xlab("Eastings") +
        ylab("Northings") +
        xlim(xlim.survey.utm) +
        ylim(ylim.survey.utm) +
        labs(fill = guide.title) +
        ggtitle(main) +
        theme_simple() +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        facet_wrap( ~ year, nrow = 2, as.table = FALSE)
    if(!is.null(guide.lims))
        g <- g + scale_fill_gradientn(colors = viridis(256),
                                      limits = guide.lims)
    else
        g <- g + scale_fill_gradientn(colors = viridis(256))

    print(g)
}



## surf_yr -------------------------------------------------
surf_yr <- function(data, method, grid,
                    exclude = TRUE,
                    exclude.method = "box",
                    exclude.dist = NULL,
                    vgm.params = NULL,
                    vgm.model = "Exp",
                    var.x = "easting_m",
                    var.y = "northing_m",
                    var.z = "temp_100",
                    var.year = "year", ...) {
    ## Fit spatial surface to data by year
    ##
    ## data = data.frame
    ## method = string specifying suface fitting method
    ## grid = data.frame with x,y columns specifying prediction grid
    ## exclude = logical, should grid cells too far from data be excluded
    ## exclude.method = see exclude_grid()
    ## exclude.dist = see exclude_grid()
    ## vgm.params = estimated parameters output from vgram_summary_yr()
    ##              only used if method == "krige"
    ## vgm.model = string specifying a gstat covariance model
    ##              only used if method == "krige"
    ## var.x = name of x-coordinate column in `data`
    ## var.y = name of y-coordinate column in `data`
    ## var.z = name of data column in `data`
    ## var.year = name of year column in `data`
    ## ... = passed to fitting method function

    years  <- unique(data[[var.year]])
    lst.yr <- vector("list", length(years))
    for(t in seq_along(years)) {
        yr <- years[t]

        tdat.t  <- data[data[[var.year]] == yr, ]
        idat.t <- data.frame(x = tdat.t[[var.x]],
                             y = tdat.t[[var.y]],
                             z = tdat.t[[var.z]])

        if(exclude) {
            excl <- exclude_grid(grid$x,
                                 grid$y,
                                 idat.t$x,
                                 idat.t$y,
                                 method = exclude.method,
                                 dist = exclude.dist)
            grd.t <- grid[!excl, ]
        } else {
            grd.t <- grid
        }

        ## Handle duplicate spatial locations
        dup.df <- data.frame(idat.t$x, idat.t$y)
        dup    <- duplicated(dup.df)
        idat.t <- idat.t[!dup, ]

        if(method == "poly")
            ss <- surf_poly(idat.t, grd.t, ...)

        if(method == "loess")
            ss <- surf_loess(idat.t, grd.t, ...)

        if(method == "idw")
            ss <- surf_idw(idat.t, grd.t, ...)

        if(method == "krige") {
            if(is.null(vgm.params)) {
                stop("vgm.params not supplied")
            }
            vg.par.t <- vgm.params[vgm.params$year == yr, ]
            vgm.t    <- gstat::vgm(psill  = vg.par.t$sill,
                                   range  = vg.par.t$range,
                                   nugget = vg.par.t$nugget,
                                   kappa  = vg.par.t$smoothness,
                                   model  = vgm.model)
            ss <- surf_krige(idat.t, grd.t, vgm = vgm.t, ...)
        }

        ss$year <- yr
        ss$grid_id <- 1:nrow(ss)
        lst.yr[[t]] <- ss

    }
    fits <- plyr::rbind.fill(lst.yr)

    return(fits)
}



## surf_cv--------------------------------------------------
surf_cv <- function(data, method,
                    vgm.params = NULL,
                    vgm.model = "Exp",
                    var.x = "easting_m",
                    var.y = "northing_m",
                    var.z = "temp_100",
                    var.year = "year", ...) {
    ## Cross validation for temp_100 surfaces
    ##
    ## This function performs a leave-one-out cross validation procedure for the
    ## specified surfacing fitting method. The function outputs a data.frame
    ## where column `z` is the "true" data and `predicted` is the LOO-CV
    ## prediciton for that data point. When method == "krige", the variogram is
    ## not recalculated for each CV iteration.
    ##
    ## data = data.frame of temp_100 data
    ## method = string, surface fitting method
    ## vgm.params = estimated parameters output from vgram_summary_yr()
    ##              only used if method == "krige"
    ## vgm.model = string specifying a gstat covariance model
    ##              only used if method == "krige"
    ## var.x = name of x-coordinate column in `data`
    ## var.y = name of y-coordinate column in `data`
    ## var.z = name of data column in `data`
    ## var.year = name of year column in `data`
    ## ... = passed to fitting method function

    years  <- unique(data[[var.year]])
    yr.lst <- vector("list", length(years))
    for(t in seq_along(years)) {
        yr <- years[t]

        ## Subset data for a single year
        tdat.t  <- data[data[[var.year]] == yr, ]
        idat.t <- data.frame(x = tdat.t[[var.x]],
                             y = tdat.t[[var.y]],
                             z = tdat.t[[var.z]])

        ## Handle duplicate spatial locations
        dup.df <- data.frame(idat.t$x, idat.t$y)
        dup    <- duplicated(dup.df)
        idat.t <- idat.t[!dup, ]

        if(method == "krige") {
                if(is.null(vgm.params)) {
                    stop("vgm.params not supplied")
                }
                vg.par.t <- vgm.params[vgm.params$year == yr, ]
                vgm.t    <- gstat::vgm(psill  = vg.par.t$sill,
                                       range  = vg.par.t$range,
                                       nugget = vg.par.t$nugget,
                                       kappa  = vg.par.t$smoothness,
                                       model  = vgm.model)
        }

        ## Loop over input data rows
        cv.lst <- vector("list", nrow(idat.t))
        for(i in 1:nrow(idat.t)) {
            dat.fit <- idat.t[-c(i), ]
            dat.rm  <- idat.t[i, ]
            dat.pr  <- dat.rm[ , c("x", "y")]

            if(method == "poly")
                ss <- surf_poly(dat.fit, dat.pr, ...)

            if(method == "loess")
                ss <- surf_loess(dat.fit, dat.pr, ...)

            if(method == "idw")
                ss <- surf_idw(dat.fit, dat.pr, ...)

            if(method == "krige")
                ss <- surf_krige(dat.fit, dat.pr, vgm = vgm.t, ...)

            ## Clean-up
            names(ss)[names(ss) == "z"] <- "predicted"
            names(ss)[names(ss) == "se.fit"] <- "se.predicted"
            ss$z <- dat.rm$z
            ss$row.ind <- i

            cv.lst[[i]] <- ss
        }
        cv.df <- plyr::rbind.fill(cv.lst)
        cv.df$error <- cv.df$z - cv.df$predicted
        cv.df$year  <- yr
        yr.lst[[t]] <- cv.df
    }
    cv.all <- plyr::rbind.fill(yr.lst)

    return(cv.all)
}



## vgram_plot_yr -------------------------------------------
vgram_plot_yr <- function(vgram, type = c("byyear", "mean"),
                          ymax = 1.5, main = "") {
    ## Plot results of vgram_fit_yr()
    ##
    ## vgram = output from vgram_fit_yr()
    ## type = "byyear" or "mean"; "byyear" plots the fitted variogram by year in
    ##         separate panels; "mean" plots the annual variograms and overlays
    ##         the mean variogram in a single panel
    ## ymax = y-axis maximum
    ## main = string for plot title

    type <- type[1]

    ## Plot variograms + fits by year
    if(type == "byyear") {
        vg.plot <- vector("list", length(vgram))
        for(i in seq_along(vgram)) vg.plot[[i]] <- vgram[[i]]$plot
        vg.plot <- plyr::rbind.fill(vg.plot)
        g <- ggplot(vg.plot) +
            geom_point(aes(x = dist, y = gamma), color = "grey50", size = 0.5) +
            geom_line(aes(x = dist, y = matern), size = 0.8, color = "red3") +
            ylab("Semivariance") +
            xlab("Distance (m)") +
            ggtitle(main) +
            theme_simple() +
            facet_wrap( ~ year, as.table = FALSE, scales = "free_y")
        print(g)
    }

    ## Plot annual fitted variograms + overlay mean
    if(type == "mean") {
        vg.yr <- vgram_summary_yr(vgram)
        vg.x <- vgram[[1]]$vgram$centers
        plot(0, 0, type = "n",
             xlim = c(0, max(vg.x)),
             ylim = c(0, ymax),
             ylab = "Semivariance",
             xlab = "Distance (m)",
             main = main,
             las = 1,
             axes = FALSE)
        box(col = "grey50")
        axis(side = 1, lwd = 0, lwd.tick = 1, col = "grey50")
        axis(side = 2, lwd = 0, lwd.tick = 1, las = 1, col = "grey50")
        for(i in 1:nrow(vg.yr)) {
            sill.i <- vg.yr$sill[i]
            rho.i  <- vg.yr$range[i]
            smoothness.i <- vg.yr$smoothness[i]
            nugget.i <- vg.yr$nugget[i]
            m0 <- matern(0, rho.i, smoothness.i, sill.i,
                         nugget = nugget.i, type = "semivariance")
            mfit <- matern(vg.x, rho.i, smoothness.i, sill.i,
                           nugget = nugget.i, type = "semivariance")
            lines(c(m0, vg.x), c(m0, mfit), col = "grey60")
        }
        sill.mu <- mean(vg.yr$sill)
        rho.mu  <- mean(vg.yr$range)
        m0 <- matern(0, rho.mu, smoothness.i, sill.mu,
                     nugget = nugget.i, type = "semivariance")
        mfit <- matern(vg.x, rho.mu, smoothness.i, sill.mu,
                       nugget = nugget.i, type = "semivariance")
        lines(c(m0, vg.x), c(m0, mfit), col = "black", lwd = 4)
    }

}

if(FALSE) {

    dat <- temp.100[temp.100$month > 5, ]
    vg   <- vgram_fit_yr(dat, years.survey)
    vgram_plot_yr(vg, type = "byyear")
    vgram_plot_yr(vg, type = "mean")

}



## vgram_summary_yr ----------------------------------------
vgram_summary_yr <- function(vgram) {
    ## Summarize vgram_fit_yr() results by year
    ##
    ## Returns a data.frame with each row summarizing a vgram fit for a
    ## particular year
    ##
    ## vgram = list ouput from vgram_fit_yr()

    vg.yr <- plyr::ldply(vgram, .id = "year", summarize,
                         range      = fit$par["range"],
                         sill       = fit$par["sill"],
                         nugget     = fit$nugget,
                         smoothness = fit$smoothness,
                         converge   = fit$convergence,
                         obj        = fit$value,
                         degree     =  fit$degree)
    return(vg.yr)
}

if(FALSE) {

    dat <- temp.100[temp.100$month > 5, ]
    vg   <- vgram_fit_yr(dat, years.survey)
    vgram_summary_yr(vg)

}



## vgram_fit_yr --------------------------------------------
vgram_fit_yr <- function(data,
                         degree = 2,
                         dmax = 200e3,
                         breaks = seq(0, dmax, 3e3),
                         smoothness = 0.5,
                         var.x = "easting_m",
                         var.y = "northing_m",
                         var.z = "temp_100",
                         var.year = "year") {
    ## Fit variogram to data by year
    ##
    ## This function uses the fields::vgram function to calculate the empirical
    ## variogram and fits a Matern covariance function to the empirical
    ## variogram. The vgram and model are calculated and fitted for each year
    ## separately.
    ##
    ## The function ouputs a list with one element for each year. Each year
    ## element is composed of a list with three elements:
    ##   $vgram = output of fields::vgram function
    ##   $fit   = output of fit_vgram() function
    ##   $plot  = a data.frame for easy plotting of the emprirical vgram and
    ##            fitted model
    ##
    ## data = data.frame, should have a "year" column
    ## degree = degree of polynomial to calc resids from
    ## dmax = maximum distance to consider; see fields::vgram
    ## breaks = breaks points of variogram bins; see fields::vgram
    ## smoothness = matern smoothness parameter
    ## var.x = name of x-coordinate column in `data`
    ## var.y = name of y-coordinate column in `data`
    ## var.z = name of data column in `data`
    ## var.year = name of year column in `data`

    deg   <- degree
    brks  <- breaks
    years <- unique(data[[var.year]])

    ## Fit variogram by year
    vg.lst <- vector("list", length(years))
    pt.lst <- vector("list", length(years))
    names(vg.lst) <- years
    for(t in seq_along(years)) {

        # browser()
        # cat(t, "\n")

        ## Prep data
        yr     <- years[t]
        dat.t  <- data[data[[var.year]] == yr, ]
        idat.t <- data.frame(x = dat.t[[var.x]],
                             y = dat.t[[var.y]],
                             z = dat.t[[var.z]])

        ## Fit variograms
        sls <- spatial::surf.ls(deg, idat.t)
        mat <- matrix(c(idat.t$x, idat.t$y), ncol = 2)
        vg  <- fields::vgram(mat, y = resid(sls), dmax = dmax, breaks = brks)
        fv  <- fit_vgram(vg$centers, vg$stats["mean", ],
                         smoothness = smoothness,
                         plot = FALSE, plot.main = yr)
        fv$degree <- deg

        ## Create df for plotting
        rho.t  <- fv$par["range"]
        sill.t <- fv$par["sill"]
        mfit   <- matern(c(0, vg$centers), rho.t, smoothness, sill.t,
                         nugget = 0, type = "semivariance")
        plt <- data.frame(year   = yr,
                          dist   = c(0, vg$centers),
                          gamma  = c(0, vg$stats["mean", ]),
                          matern = mfit,
                          range  = rho.t,
                          sill   = sill.t,
                          nugget = 0,
                          degree = deg)

        ## Prepare outputs
        vg.lst[[t]] <- list(vgram = vg, fit = fv, plot = plt)

        rm(vg); gc()

    }
    return(vg.lst)
}

if(FALSE) {

    dat <- temp.100[temp.100$month > 5, ]

    vg <- vgram_fit_yr(dat, years.survey)
    names(vg)
    vg[[1]]$fit

}



## surf_*() ------------------------------------------------
##
## These functions fit surfaces to spatial data:
##
## surf_poly:   fit global polynomial surface
## surf_loess:  fit local polynomial surface
## surf_idw:    interpolate surface using inverse distance weighting
## surf_krige:  interpolate surface using kriging
##
## Each function returns a data.frame with columns:
##   $x       x-coordinates of prediction surface
##   $y       y-coordinates of prediction surface
##   $z       value of variable at location x,y
##   $label   unique identifier for surface fitting method
##
##
## Common arguments across functions:
##    data = data.frame with columns x,y,z
##    grid = data.frame with columns x,y, is used as the prediction grid for
##           the smooth surface
##    label = optional label for method, added to output df as a column
##    degree = degree of polynomial
##    power = power for inverse distance weights
##    se = logical, should standard error be calculated
##    span = see ?loess
##     ... = passed to model fitting function

surf_poly <- function(data, grid, degree = 2, label = NULL, ...) {
    fit    <- spatial::surf.ls(degree, data, ...)
    pr.mat <- predict(fit, x = grid$x, y = grid$y)
    pr.df  <- data.frame(x = grid$x, y = grid$y, z = pr.mat)
    pr.df$method <- "poly"
    pr.df$degree <- degree
    pr.df$label <- label
    if(is.null(label)) {
        pr.df$label <- paste0("poly-deg", degree)
    }
    return(pr.df)
}


surf_loess <- function(data, grid, se = FALSE,
                       degree = 2, span = 1, label = NULL, ...) {
    fit.loess <- loess(z ~ x * y, data, degree = degree, span = span,
                       control = loess.control(surface = "direct"), ...)
    pr.loess <- predict(fit.loess, grid, se = se)
    pr.df <- data.frame(x = grid$x, y = grid$y)
    if(se) {
        pr.df$z <- pr.loess$fit
        pr.df$se.fit <- pr.loess$se.fit
    } else {
        pr.df$z <- pr.loess
    }
    pr.df$method <- "loess"
    pr.df$degree <- degree
    pr.df$span   <- span
    pr.df$label  <- label
    if(is.null(label)) {
        pr.df$label <- paste0("loess-deg", degree, "-span", span)
    }
    return(pr.df)
}


surf_idw <- function(data, grid, power = 2, label = NULL,
                     debug.level = 0, ...) {
    fit.idw <- gstat::idw(z ~ 1, data = data,
                          locations = ~ x + y,
                          newdata = grid, idp = power,
                          debug.level = debug.level, ...)
    pr.idw <- data.frame(x = fit.idw$x, y = fit.idw$y, z = fit.idw$var1.pred)
    pr.idw$method <- "idw"
    pr.idw$degree <- power
    pr.idw$label <- label
    if(is.null(label)) {
        pr.idw$label <- paste0("idw-power", power)
    }
    return(pr.idw)
}

surf_krige <- function(data, grid, vgm, formula = NULL, degree = 0,
                       label = NULL, debug.level = 0, ...) {
    ## vgm = gstat::vgm output
    ##
    ## if formula is non-null, degree is ignored, otherwise, degree specifies
    ## the degree of polynomial and should be in c(0, 1, 2).

    if(is.null(formula)) {
        if(degree == 0) fo <- z ~ 1
        if(degree == 1) fo <- z ~ x + y
        if(degree == 2) fo <- z ~ x + y + I(x*x) + I(y*y) + I(x*y)
    } else {
        fo <- formula
    }
    fit.kr <- gstat::krige(fo, locations = ~ x + y,
                           data = data, model = vgm, newdata = grid,
                           debug.level = debug.level, ...)
    pr.kr <- data.frame(x = fit.kr$x, y = fit.kr$y, z = fit.kr$var1.pred)
    pr.kr$se.fit <- sqrt(ifelse(fit.kr$var1.var < 0, 0, fit.kr$var1.var))
    pr.kr$method <- "krige"
    pr.kr$degree <- degree
    pr.kr$label <- label
    if(is.null(label)) {
        if(is.null(formula))
            pr.kr$label <- paste0("krige-deg", degree)
        else
            pr.kr$label <- paste0("krige")
    }
    return(pr.kr)
}

if(FALSE) {

    dat <- temp.100[temp.100$month > 5, ]
    dat.2015 <- dat[dat$year == 2015, ]
    idat.2015 <- data.frame(x = dat.2015$easting_m,
                            y = dat.2015$northing_m,
                            z = dat.2015$temp_100)

    ss <- surf_poly(idat.2015, grd.excl)

    ss <- surf_loess(idat.2015, grd.excl, se = FALSE)
    ss <- surf_loess(idat.2015, grd.excl, se = TRUE)

    ss <- surf_idw(idat.2015, grd.excl, power = 2)

    fo  <- z ~ x + y + I(x*x) + I(y*y) + I(x*y)
    fv  <- vgm(psill = 0.7, model = "Exp", range = 8199, nugget = 0)
    ss1 <- surf_krige(idat.2015, grd.excl, vgm = fv, degree = 2)
    ss2 <- surf_krige(idat.2015, grd.excl, vgm = fv, formula = fo)
    ss2 <- surf_krige(idat.2015, grd.excl, vgm = fv, formula = fo, maxdist = 150e3)

}



## exclude_grid --------------------------------------------
exclude_grid <- function(x.grid, y.grid, x.data, y.data,
                         method = "box",
                         dist = (x.grid[2] - x.grid[1]) / 2,
                         plot = FALSE, ...) {
    ## Exclude grid cells that are too far from data
    ##
    ## Given x,y coordinate data, grid cells will be excluded if they are futher
    ## than `dist` from an x,y data point. For large datasets, consider using
    ## `method = "unit"` to speed up the exclusion algorithm.
    ##
    ## exclude.methods:
    ##   euclidean = only include grid cells that have a data point within
    ##               dist. Euclidean distances are used. dist
    ##               represents the radius of a circle around the grid cell
    ##               center in which to search for data points
    ##   box = only include grid cells that have data +/- dist from the
    ##         grid cell center. dist is the width of a box around the
    ##         center in which to search for data points
    ##   unit = see mgcv::exclude.too.far
    ##   lonlat = like euclidean, but exclude grid cells based on dist (in
    ##            meters) when input is lon,lat coordinates. The haversine
    ##            equation is used to calc distances
    ##
    ## x.grid = x coordinate for grid cell centers
    ## y.grid = y coordinate for grid cell centers
    ## x.data = x coordinate for data
    ## y.data = y coordinate for data
    ## method = c("euclidean", "unit", "box")
    ## dist = see above
    ## plot = logical, should grid be plotted, light grey box gives bounding box
    ##        for x,y
    ## ... = passed to plot()

    res <- x.grid[2] - x.grid[1]
    grd <- data.frame(x = x.grid, y = y.grid)

    if(method == "euclidean") {
        excl <- rep(NA, nrow(grd))
        for(i in 1:nrow(grd)) {
            x.sq <- (grd[i, "x"] - x.data)^2
            y.sq <- (grd[i, "y"] - y.data)^2
            dst  <- sqrt(x.sq + y.sq)
            excl[i] <- !any(dst < dist)
        }
    }

    if(method == "box") {
        excl <- rep(NA, nrow(grd))
        for(i in 1:nrow(grd)) {
            xi.max <- grd[i, "x"] + dist
            xi.min <- grd[i, "x"] - dist
            yi.max <- grd[i, "y"] + dist
            yi.min <- grd[i, "y"] - dist
            excl[i] <- !any(x.data <= xi.max & x.data >= xi.min &
                            y.data <= yi.max & y.data >= yi.min)
        }
    }

    if(method == "lonlat") {
        excl <- rep(NA, nrow(grd))
        ll.dat <- matrix(c(x.data, y.data), ncol = 2)
        for(i in 1:nrow(grd)) {
            ll.grd <- c(grd$x[i], grd$y[i])
            dist.i <- geosphere::distHaversine(ll.grd, ll.dat)
            min.dist <- dist.i[which.min(dist.i)]
            excl[i] <- min.dist > dist
        }
    }

    if(method == "unit") {
        excl <- mgcv::exclude.too.far(grd$x, grd$y, x.data, y.data,
                                      dist = dist)
    }

    grd.excl <- grd[!excl, ]
    row.names(grd.excl) <- NULL

    if(plot) {
        plot(grd.excl$x, grd.excl$y, cex = 0.2, col = "tomato",
             pch = 21, asp = 1, ...)
        for(i in 1:nrow(grd.excl)) {
            rect(xleft = grd.excl$x[i] - (res / 2),
                 ybottom = grd.excl$y[i] - (res / 2),
                 xright = grd.excl$x[i] + (res / 2),
                 ytop = grd.excl$y[i] + (res / 2),
                 border = "grey60", lty = 1)
        }
    }

    return(excl)
}

if(FALSE) {

    x <- runif(100, min = 0, max = 100)
    y <- runif(100, min = 0, max = 100)
    plot(x, y)

    grd <- pred_grid(x, y, res = 5, plot = TRUE)
    points(x, y, cex = 0.5, col = 1, pch = 1)
    grd.excl <- exclude_grid(grd$x, grd$y, x, y, plot = TRUE)
    points(x, y, cex = 0.8, col = 1, pch = 1)

}



## pred_grid -----------------------------------------------
pred_grid <- function(x, y, res,
                      pad = 0,
                      plot = FALSE) {
    ## Create prediction grid from x,y data
    ##
    ## Given x,y coordinate data this function creates a prediction grid based
    ## on the bounding box of the x,y data. The function returns a data.frame
    ## with columns x,y of the grid cell centers. If pad=0 then all grid cell
    ## centers are located within the bounding box of x,y. The grid is created
    ## with cell dims (res x res) where res is in the same units as x,y.
    ##
    ## x = x coordinate
    ## y = y coordinate
    ## res = grid resolution in same units as x,y (res x res grid cells)
    ## pad = pad bounding box of data for grid
    ## plot = logical, should grid be plotted, light grey box gives bounding box
    ##        for x,y

    if(missing(res))
        stop("`res` not set")

    xmin <- min(x) - pad
    xmax <- max(x) + pad
    ymin <- min(y) - pad
    ymax <- max(y) + pad

    ## Create initial seq, this sequence has the undesired effect of not being
    ## centered in the bounding box of the data b/c the sequence always starts
    ## at the minimum and has some remainder on the max end of the seq
    x.grd1 <- seq(xmin, xmax, res)
    y.grd1 <- seq(ymin, ymax, res)

    ## Find remainder distance from seq calculation
    xrem <- xmax - max(x.grd1)
    yrem <- ymax - max(y.grd1)

    ## Add 1/2 of remainder to min values to center grid
    ## in bounding box
    xmin <- xmin + (xrem / 2)
    ymin <- ymin + (yrem / 2)

    ## Create centered grid sequence
    x.grd <- seq(xmin, xmax, res)
    y.grd <- seq(ymin, ymax, res)

    grd <- expand.grid(x = x.grd, y = y.grd, KEEP.OUT.ATTRS = FALSE)

    if(plot) {
        plot(grd$x, grd$y, cex = 0.2, col = "tomato", pch = 21, asp = 1)
        rect(xleft = min(x), ybottom = min(y),
             xright = max(x), ytop = max(y),
             border = "grey30", lty = 1)
        for(i in 1:nrow(grd)) {
            rect(xleft = grd$x[i] - (res / 2),
                 ybottom = grd$y[i] - (res / 2),
                 xright = grd$x[i] + (res / 2),
                 ytop = grd$y[i] + (res / 2),
                 border = "grey80", lty = 1)
        }
    }

    return(grd)
}

if(FALSE) {

    x <- runif(100, min = 0, max = 100)
    y <- runif(100, min = 0, max = 100)
    plot(x, y)

    grd <- pred_grid(x, y, res = 5, plot = TRUE)
    points(x, y, cex = 0.5, col = 1, pch = 1)

    grd <- pred_grid(x, y, res = 1, plot = TRUE)
    points(x, y, cex = 0.5, col = 1, pch = 1)

    grd <- pred_grid(x, y, plot = TRUE)
    points(x, y, cex = 0.5, col = 1, pch = 1)
}



## fit_vgram -----------------------------------------------
fit_vgram <- function(x, y,
                      smoothness = 0.5,
                      init.range = NULL,
                      init.sill = NULL,
                      nugget = 0,
                      plot = FALSE,
                      plot.main = "",
                      xlim = NULL,
                      ylim = NULL,
                      optim.control = list(), ...) {
    ## Fit a matern covariance function to a semi-variogram
    ##
    ## This function fits a matern covariance function to semi-variogram. The
    ## matern covariance function has two parameters that are fitted: range and
    ## sill. Nugget component is not implemented. A warning is issued if optim()
    ## fails to converge in any manner. The function returns the output from
    ## optim().
    ##
    ## x = vector of variogram distances
    ## y = vector of semivariances
    ## init.range = initial value for range parameter
    ## init.sill = initial value for sill parameter
    ## nugget = fixed value for nugget
    ## plot = logical, should fitted model be plotted
    ## plot.main = string to use for plot title
    ## optim.control = list of controls passed to optim; see ?optim
    ## ... = passed to optim()

    if(is.null(init.range)) {
        init.range <- max(x) / 3
    }

    if(is.null(init.sill)) {
        init.sill <- mean(y[(length(y)-5):length(y)])
    }

    obj <- function(params, x, y, smoothness, nugget) {
        if(any(c(params["range"] <= 0, params["sill"] <= 0))) {
            ssq <- Inf
        } else {
            co <- matern(x, smoothness = smoothness,
                         range = params["range"],
                         sill = params["sill"],
                         nugget = nugget,
                         type = "semivariance")
            resid <- y - co
            ssq   <- sum(resid * resid)  ## OLS
        }
        return(ssq)
    }

    inits <- c(range = init.range, sill = init.sill)
    op <- optim(inits, obj, x = x , y = y,
                smoothness = smoothness, nugget = nugget,
                control = optim.control, ...)

    if(plot) {
        sill <- op$par["sill"]
        rho  <- op$par["range"]
        m0   <- matern(0, rho, smoothness, sill, nugget, "semivariance")
        mfit <- matern(x, rho, smoothness, sill, nugget, "semivariance")
        plot(x, y,
             main = plot.main,
             xlim = xlim,
             ylim = ylim,
             ylab = paste("semivariance, sill =", round(sill, 2)),
             xlab = paste("distance, range =", round(rho, 2)))
        lines(c(m0, x), c(m0, mfit), col = 2)
    }

    if(op$convergence != 0) {
        warning("Optim did not converge")
    }

    op$inits      <- inits
    op$nugget     <- nugget
    op$smoothness <- smoothness

    return(op)
}

if(FALSE) {

    library(spatial)
    data(topo, package = "MASS")
    topo.kr <- surf.ls(2, topo)
    vg <- spatial::variogram(topo.kr, nint = 25)

    fit_vgram(vg$x, vg$y, smoothness = 0.5, plot = TRUE)
    fit_vgram(vg$x, vg$y, smoothness = 1.0, plot = TRUE)
    fit_vgram(vg$x, vg$y, smoothness = 1.5, plot = TRUE)

}



## matern --------------------------------------------------
matern <- function(distance,
                   range = 1,
                   smoothness = 1,
                   sill = 1,
                   nugget = 0,
                   type = c("covariance", "semivariance")) {
    ## Matern covariance function (isotropic)
    ##
    ## This function evaluates the matern covariance function over the supplied
    ## distances. The function returns a vector of matern covariance values. If
    ## type is "covariance" a matern correlation function is used, if type is
    ## "semivariance" the correlation function is transformed to semivariances,
    ## possibly with a nugget. The nugget argument is ignored when type is
    ## "covariance".
    ##
    ## distance = vector of distances between points
    ## range = (rho or kappa) decorrelation rate
    ##  - controls the rate of decay with distance, with larger values of rho
    ##    corresponding to more highly correlated observations
    ## smoothness = (nu) VAST fixes this at 1, when smoothness = 0.5 the
    ##              function corresponds to an exponential covariance
    ## sill = marginal variance
    ## nugget = y-intercept, only used if type = "semivariance"
    ##
    ## Modified from fields::Matern

    d <- distance
    if (any(d < 0))
        stop("distance argument must be nonnegative")
    d[d == 0] <- 1e-10

    type <- type[1]

    num <- (d / range)^smoothness
    den <- gamma(smoothness) * (2^(smoothness - 1))
    cov <- (num / den) * besselK((d / range), smoothness)
    if(type == "covariance") {
        ret <- sill * cov
    }
    if(type == "semivariance") {
        ret <- (sill - nugget) * (1 - cov) + nugget
    }
    return(ret)
}

if(FALSE) {

    d <- seq(0, 5, 0.1)
    m <- matern(d, sill = 1, nugget = 0)
    plot(d, m, ylim = c(0, 1))

    m <- matern(d, sill = 1, nugget = 0, type = "semivariance")
    plot(d, m, ylim = c(0, 1))

    library(spatial)
    data(topo, package = "MASS")

    topo.kr <- surf.ls(2, topo)
    d <- seq(0, 7, 0.1)
    corr <- spatial::correlogram(topo.kr, 25)
    lines(d, matern(d, range = 0.5, smoothness = 1, sill = 1), col = 2)

    corr <- spatial::variogram(topo.kr, 25)
    lines(d, matern(d, range = 0.7, smoothness = 1, sill = 1e3, type = "semivariance"), col = 2)

    ## Change marginal variance
    corr <- correlogram(topo.kr, 25)
    lines(d, matern(d, range = 0.5, smoothness = 1, sill = 0.5, nugget = 0), col = 2)

    ## Exponential form (smoothness = 0.5)
    corr <- correlogram(topo.kr, 25)
    lines(d, matern(d, range = 0.5, smoothness = 0.5), col = 2)
    lines(d, expcov(d, d = 0.5), col = 3, lty = 2)

}





## check_profiles ------------------------------------------
check_profiles <- function(profiles, events, verbose = FALSE) {
    ## Sanity data checks for profiles+events data.frames's
    ##
    ## Function returns TRUE if no errors were found
    ##
    ## profiles = data.frame of profiles
    ## events = data.frame of events (one row per profile)
    ## verbose = print potentially useful info

    if(verbose) cat("Checking for multiples...", "\n")
    tst <- is.null(temp_find_multiples(events,
                                       select = c("lat", "lon", "date")))
    if(!isTRUE(tst)) {
        stop("Multiples found in `events`")
    }

    if(verbose) cat("Checking for duplicate keys...", "\n")
    tst <- sum(duplicated(events$profile_key)) == 0
    if(!isTRUE(tst)) {
        stop("Duplicate `profile_key` in `events`")
    }

    if(verbose) cat("Comparing keys in events and profiles...", "\n")
    tst <- all.equal(order(events$profile_key),
                     order(unique(profiles$profile_key)))
    if(!isTRUE(tst)) {
        stop("Error in `profile_key` between `events` and `profiles`")
    }

    if(verbose) cat("Done!", "\n")
    return(TRUE)
}

## add_profile_category ------------------------------------
add_profile_category <- function(data, depth,
                                 max.depth.var = "max_depth_m",
                                 bottom.depth.var = "bottom_depth_m",
                                 category.var = "category_profile") {
    ## Assign a category to water column profiles
    ##
    ## This function takes as input a data.frame of water column profiles and
    ## assigns each profile to one of four categories:
    ##   Category 1: max profile depth > `depth` & bottom depth > `depth`
    ##   Category 2: max profile depth < `depth` & bottom depth < `depth`
    ##   Category 3: max profile depth < `depth` & bottom depth > `depth`
    ##   Category 4: max profile depth > `depth` & bottom depth < `depth`
    ##
    ## The input data.frame should have columns specifying the max depth of the
    ## profile and the bottom depth at the location of the profile. The function
    ## outputs the input data.frame with an additional column `category.var`
    ##
    ## data = a data.frame with profiles to categorize
    ## depth = threshold depth for categorizing
    ## categories = optional vector of length four specifying names for each of
    ##              the four categories defined above.
    ## max.depth.var = column name in `data` giving the max depth of the profile
    ## bottom.depth.var = column name in `data` giving the bottom depth
    ##                    for profile
    ## category.var = name of the category column to append in the output

    categories <- categories_define(depth)

    dat <- data.table::copy(data)
    dat[[category.var]] <- NA

    dat[[category.var]] <- ifelse(dat[[max.depth.var]] >= depth &
                                         dat[[bottom.depth.var]] >= depth,
                                         categories[1],
                                         dat[[category.var]])
    dat[[category.var]] <- ifelse(dat[[max.depth.var]] < depth &
                                         dat[[bottom.depth.var]] < depth,
                                         categories[2],
                                         dat[[category.var]])
    dat[[category.var]] <- ifelse(dat[[max.depth.var]] < depth &
                                         dat[[bottom.depth.var]] >= depth,
                                         categories[3],
                                         dat[[category.var]])
    dat[[category.var]] <- ifelse(dat[[max.depth.var]] >= depth &
                                         dat[[bottom.depth.var]] < depth,
                                         categories[4],
                                         dat[[category.var]])
    return(dat)
}



## categories_define ---------------------------------------
categories_define <- function(depth) {
    ## Create vector of profile categories based on depth
    ##
    ## depth = depth for category vector

    cat1 <- paste0("profile > ", depth, "m & bottom > ", depth, "m")
    cat2 <- paste0("profile < ", depth, "m & bottom < ", depth, "m")
    cat3 <- paste0("profile < ", depth, "m & bottom > ", depth, "m")
    cat4 <- paste0("profile > ", depth, "m & bottom < ", depth, "m")
    categories <- c(cat1, cat2, cat3, cat4)
    return(categories)
}

if(FALSE) {
    categories_define(10)
}



## temp_find_multiples -------------------------------------
temp_find_multiples <- function(events, select) {
    ## Find multiple temperature profiles
    ##
    ## This function returns a data.frame of the rows in events that have
    ## non-unique combination of the variables specified by `select`. A new
    ## column, `multiple_id` is also created which provides a unique ID for each
    ## mulitple combination. E.g., if select = c("lon", "lat", "date") then any
    ## profiles (i.e., rows) in events that occur at the identical location and
    ## date are returned with all non-unique combinations row binded into a
    ## data.frame.
    ##
    ## events = summary of profiles data.frame w/ one row per profile
    ## select = columns names to use in finding multiples

    ## Find Mulitples
    dups <- plyr::ddply(events, .(year), function(x) {
                        ll <- subset(x, select = select)
                        dup.ind <- which(duplicated(ll) |
                                         duplicated(ll, fromLast = TRUE))
                        dup <- x[dup.ind, ]
                        out <- dup[order(dup$lon, dup$lat), ]
                        out
                    })

    ## Add a unique ID for each duplicate combination
    gr.sp <- plyr::dlply(dups, select)
    gr.ls <- lapply(seq_along(gr.sp), function(x) {
               dat.event <- gr.sp[[x]]
               dat.event$multiple_id <- x
               dat.event
           })
    gr <- plyr::rbind.fill(gr.ls)

    return(gr)
}



## ocean_hake_process --------------------------------------
ocean_hake_process <- function(data, years = NULL, source = NULL) {
    ## Process hake survey oceanographic data
    ##
    ## data = data.frame of water column profiles (in long format)
    ## years = if non-null, vector of years to subset
    ## source = if non-null, add a $source column set to this value

    dat <- data.table::copy(data)

    ## Subset years
    if(!is.null(years)) {
        dat <- dat[dat$year %in% years, ]
    }

    ## Set source
    if(!is.null(source)) {
        dat$source <- source
    }

    ## Specify variants of column names
    depth.cols <- c("depth", "Depth_m", "depth_m", "Depth", "DepSM")
    temp.cols <- c("temperature", "Temp_C", "Temperature_degC",
                   "temperature_degC", "Temperatur", "Temp1", "temperature_c")
    salininity.cols <- c("salinity", "Salinity_PSU", "Salinity_psu")
    oxygen.cols1 <- c("Oxygen_ml_per_l", "DO_mL_L")
    oxygen.cols2 <- c("Oxygen_ml_per_l_Secondary")

    ## Specify columns to subset
    common.cols <- c("year",
                     "ship",
                     "survey",
                     "instrument",
                     "profile_id",
                     "source",
                     "country",
                     "date",
                     "time",
                     "lon",
                     "lat",
                     "bottom_depth_m",
                     "dataset",
                     depth.cols,
                     temp.cols,
                     salininity.cols,
                     oxygen.cols1,
                     oxygen.cols2,
                     "n_bin",
                     "u_velocity_m_s",
                     "v_velocity_m_s")
    nam <- names(dat)[names(dat) %in% common.cols]
    if(is(data, "data.table")) {
        dat.sub <- dat[ , ..nam]
    } else {
        dat.sub <- dat[ , nam]
    }

    ## Rename depth, temp, and salinity columns
    names(dat.sub)[names(dat.sub) %in% depth.cols] <- "depth_m"
    names(dat.sub)[names(dat.sub) %in% temp.cols] <- "temp_c"
    names(dat.sub)[names(dat.sub) %in% salininity.cols] <- "salinity_psu"
    names(dat.sub)[names(dat.sub) %in% oxygen.cols1] <- "oxygen_ml_l"
    names(dat.sub)[names(dat.sub) %in% oxygen.cols2] <- "oxygen2_ml_l"
    names(dat.sub)[names(dat.sub) == "bottom_depth_m"] <- "instrument_bottom_depth_m"

    ## Parse date
    d.date <- as.Date(dat.sub$date)
    dat.sub$date <- d.date
    dat.sub$month <- as.numeric(format(d.date, "%m"))
    dat.sub$day <- as.numeric(format(d.date, "%d"))
    dat.sub$jday <- as.numeric(format(d.date, "%j"))

    ## Make sure depth is an integer
    dat.sub[["depth_m"]] <- floor(dat.sub[["depth_m"]])

    ## Clean-up
    row.names(dat.sub) <- NULL

    return(dat.sub)
}


## add_bottom_depth ----------------------------------------
add_bottom_depth <- function(data,
                             max.search,
                             min.depth,
                             etopo.data,
                             var.lon = "lon",
                             var.lat = "lat",
                             etopo.lon = "lon",
                             etopo.lat = "lat",
                             var.bdepth = "bottom_depth_m",
                             var.order = NULL,
                             var.distance = NULL) {
    ## Add bottom depth to the input data.frame
    ##
    ## This function takes as input a data.frame with lon and lat columns
    ## specified by `var.lon` and `var.lat` and (1) matches bottom depth to each
    ## locations, and (2) appends the bottom depth to the input data.frame as a
    ## new column.
    ##
    ## data = a data.frame with lon and lat columns
    ## max.search = see `get_etopo_depth`
    ## min.depth = see `get_etopo_depth`
    ## etopo.data = see `get_etopo_depth`
    ## var.lon = longitude column name in `data`
    ## var.lat = latitude column name in `data`
    ## etopo.lon = longitude column name in `etopo.data`
    ## etopo.lat = latitude column name in `etopo.data`
    ## var.bdepth = name of new column for bottom depth
    ## var.order = if non-null the output data.frame is ordered by this column
    ## var.distance = if non-null, specifies the column name to append the
    ##                distance in which the bottom depth was from the location

    if(var.bdepth %in% names(data))
        stop("`var.bdepth` already exists in data")
    if(!var.lon %in% names(data))
        stop("`var.lon` not found")
    if(!var.lat %in% names(data))
        stop("`var.lat` not found")
    if(!etopo.lon %in% names(etopo.data))
        stop("`etopo.lon` not found")
    if(!etopo.lat %in% names(etopo.data))
        stop("`etopo.lat` not found")

    dat <- data.table::copy(data)

    ## Subset data lon,lat
    loc.dat <- matrix(c(dat[[var.lon]], dat[[var.lat]]), ncol = 2)

    ## Match locations with bottom depth
    bdepth <- get_etopo_depth(loc.dat,
                              etopo.data = etopo.data,
                              max.search = max.search,
                              min.depth = min.depth)

    ## Clean-up depth df
    if(!is.null(var.distance)) {
        bdepth[[var.distance]] <- bdepth$distance_m
    }
    bdepth$distance_m <- NULL

    if(var.bdepth != "bottom_depth_m") {
        bdepth[[var.bdepth]] <- bdepth$bottom_depth_m
        bdepth$bottom_depth_m <- NULL
    }

    ## Handle duplicate entries
    dups <- which(duplicated(bdepth))
    if(length(dups) > 0) {
        bdepth <- bdepth[-dups, ]
    }

    ## Set columns to merge by
    names(bdepth)[names(bdepth) == etopo.lon] <- var.lon
    names(bdepth)[names(bdepth) == etopo.lat] <- var.lat
    m.by <- intersect(names(dat), names(bdepth))

    ## Merge data and bottom.depth
    dat.merge <- merge(dat, bdepth,
                       by = m.by,
                       all.x = TRUE,
                       sort = FALSE)

    ## Order data
    if(!is.null(var.order)) {
        dat.merge <- dat.merge[order(dat.merge[[var.order]]), ]
    }

    ## Clean-up
    row.names(dat.merge) <- NULL

    return(dat.merge)
}



## download_etopo ------------------------------------------
download_etopo <- function(write.path, verbose = TRUE) {

    ## Download ETOPO1 global relief model data
    ##
    ## This function downloads ETOPO1 data and saves it to disk in netCDF
    ## format. The data set is large and can take a while download. The
    ## downloaded data are in a gzipped netCDF format.
    ##
    ## See: https://www.ngdc.noaa.gov/mgg/global/global.html
    ##
    ## write.path = directory to save data file in, e.g., "./data/"

    if(verbose)
        cat("Downloading ETOPO1 data...", "\n")

    if(!dir.exists(write.path))
        dir.create(write.path, recursive = TRUE)

    url1      <- "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/"
    url2      <- "ice_surface/grid_registered/netcdf/"
    f.name    <- "ETOPO1_Ice_g_gmt4.grd.gz"
    b.url     <- paste0(url1, url2, f.name)
    dest.file <- paste0(write.path, f.name)
    download.file(b.url, destfile = dest.file)
}



## get_etopo_depth -----------------------------------------
get_etopo_depth <- function(locations, etopo.data,
                            max.search = 0.1,
                            min.depth = 20) {

    ## Find the ETOPO grid cell closest to a specified lon,lat location that has
    ## a depth at least as deep as a threshold depth
    ##
    ## Arguments:
    ##  locations = two column matrix of locations, first column should be
    ##              longitude, second column should be latitude
    ##  etopo.data = ETOPO data output from etopo_load()
    ##  max.search = value in which to expand the `lon` and `lat` values. Used
    ##               to determine the size of the ETOPO region (i.e., number of
    ##               grid cells) to search for the nearest grid cell center
    ##  min.depth = minimum depth to accept for the given location
    ##
    ## This function returns a data.frame with one row for each input location
    ## and four columns with the lon and lat of the grid cell center chosen, the
    ## depth at that grid cell, and the distance from that grid cell center and
    ## the input lon,lat location.
    ##
    ## Subsetting the relief matrix does two things:
    ##   1. Reduces the number of grid cells to calculate the distance for
    ##   2. Outputs a *square* matrix with equal number of lat and lon values,
    ##      having a square matrix is necessary to leverage the `outer()`
    ##      function

    elon.range <- range(etopo.data$lon)
    elat.range <- range(etopo.data$lat)

    n <- nrow(locations)
    out.mat <- matrix(NA, nrow = n, ncol = 4)
    for(i in 1:n) {

        lon <- locations[i, 1]
        lat <- locations[i, 2]

        ## Extract data
        elon  <- etopo.data$lon
        elat  <- etopo.data$lat
        etopo <- etopo.data$depth

        ## Subset lon and lat
        elat.min <- lat - max.search
        elat.max <- lat + max.search
        elon.min <- lon - max.search
        elon.max <- lon + max.search

        if(elat.min < elat.range[1] | elat.max > elat.range[2])
            stop("Expanded location lat outside of `etopo.data` geographic range")

        if(elon.min < elon.range[1] | elon.max > elon.range[2])
            stop("Expanded location lon outside of `etopo.data` geographic range")

        elon.sub.ind <- which(elon >= elon.min & elon <= elon.max)
        elat.sub.ind <- which(elat >= elat.min & elat <= elat.max)
        elon.sub <- elon[elon.sub.ind]
        elat.sub <- elat[elat.sub.ind]

        ## Subset relief matrix
        etopo.sub <- etopo[elon.sub.ind, elat.sub.ind]
        etopo.sub <- ifelse(etopo.sub < min.depth, NA, etopo.sub)

        ## Calculate distances between relief cells and input lon,lat
        dist <- outer(1:nrow(etopo.sub), 1:ncol(etopo.sub), function(r, c) {
            ll <- matrix(c(elon.sub[r], elat.sub[c]), ncol = 2)
            geosphere::distHaversine(c(lon, lat), ll)
        })

        ## Set dist to NA when depth is NA
        dist.na <- ifelse(is.na(etopo.sub), NA, dist)

        ## Find matrix indices with minimum distances
        ## This can result in more than one relief cell having
        ## the same minimum distance
        mat.ind <- which(dist.na == min(dist.na, na.rm = TRUE), arr.ind = TRUE)

        ## If there is only one minimum distance, use it, otherwise use the
        ## relief cell with the deepest depth
        if(nrow(mat.ind) == 1) {
            out.depth <- etopo.sub[mat.ind[1], mat.ind[2]]
            ## out.lat <- elat.sub[mat.ind[2]]
            ## out.lon <- elon.sub[mat.ind[1]]
            out.dist <- dist[mat.ind[1], mat.ind[2]]
        } else {
            depths <- apply(mat.ind, 1, function(x) etopo.sub[x[1], x[2]])
            out.ind <- which(depths == min(depths)) ## get deepest depth
            ind <- mat.ind[out.ind, ]
            out.depth <- etopo.sub[ind[1], ind[2]]
            ## out.lat <- elat.sub[ind[2]]
            ## out.lon <- elon.sub[ind[1]]
            out.dist <- dist[ind[1], ind[2]]
        }

        ## out.mat[i, ] <- c(out.lon, out.lat, out.depth, out.dist)
        out.mat[i, ] <- c(lon, lat, out.depth, out.dist)
    }
    colnames(out.mat) <- c("lon", "lat", "bottom_depth_m", "distance_m")
    out.df <- as.data.frame(out.mat)
    return(out.df)
}



## etopo_load ----------------------------------------------
etopo_load <- function(path,
                       min.depth = 0,
                       lat.min = 30,
                       lat.max = 60,
                       lon.min = -140,
                       lon.max = -115) {

    ## Load ETOPO1 data into R
    ##
    ## This function loads ETOPO1 data in netCDF format into R. The
    ## function outputs a list with three elements:
    ##  1. $lon    vector of longitudes
    ##  2. $lat    vector of latitudes
    ##  5. $depth  matrix of relief data (in meters)
    ##
    ## The ETOPO matrix has dimensions:
    ##      lon (rows) x lat (columns)
    ##
    ## Depths are converted so that negative values are above the sea surface,
    ## as per the convention in oceanography
    ##
    ## Lat and Lon give the center of the grid cell
    ##
    ## path = path to ETOPO1 data file
    ## min.depth = values greater than `min.depth` are set to NA
    ## lat.min, lon.min = minimum lat and lon values to retain
    ## lat.max, lon.max = maximum lat and lon values to retain

    R.utils::gunzip(path, remove = FALSE)
    path.nogz <- gsub("\\.gz$", "", path)

    lat.name   <- "y"
    lon.name   <- "x"
    etopo.name <- "z"

    etopo.ncdf <- ncdf4::nc_open(path.nogz)
    dat.lat    <- ncdf4::ncvar_get(etopo.ncdf, varid = lat.name)
    dat.lon    <- ncdf4::ncvar_get(etopo.ncdf, varid = lon.name)
    dat.etopo  <- ncdf4::ncvar_get(etopo.ncdf, varid = etopo.name)
    ncdf4::nc_close(etopo.ncdf)

    ## replace missing values
    dat.etopo[dat.etopo < -50000] <- NA

    ## convert negatives
    dat.etopo <- dat.etopo * -1

    ## min.depth
    dat.etopo[dat.etopo < min.depth] <- NA

    ## spatially subset data
    lat.sub.ind <- which(dat.lat >= lat.min & dat.lat <= lat.max)
    lon.sub.ind <- which(dat.lon >= lon.min & dat.lon <= lon.max)
    lst <- list(lon   = dat.lon[lon.sub.ind],
                lat   = dat.lat[lat.sub.ind],
                depth = dat.etopo[lon.sub.ind, lat.sub.ind])

    unlink(path.nogz)

    return(lst)
}



## dist_min ------------------------------------------------
dist_min <- function(p1, p2,
                     method = "euclidean",
                     max.dist = NULL,
                     max.days = NULL,
                     unique = FALSE, 
                     loud = FALSE) {
    ## Match points based on minimum distances
    ##
    ## This function takes as input two sets of x,y locations and for each point
    ## in the first set (p1), finds the point in the second set (p2) that
    ## is closest (i.e., has the minimum distance). The function returns a list
    ## with five elements:
    ##   $index = vector the same length as nrow(p1) that gives the index of p2
    ##            that has the minimum distance
    ##   $dist = vector the same length as nrow(p1) giving distance between p1
    ##            and closest p2
    ##   $max.dist = value of `max.dist` input
    ##   $max.days = value of `max.days` input
    ##   $unique = value of `unique` input
    ##
    ## Values in $index and $dist are set to NA if no p2 point was
    ## assigned to a particular p1 observation
    ##
    ## The `max.days` and `max.dist` arguments can be used to refine how to
    ## match the p2 data to p1 observations. Settings `max.days` to a non-null
    ## value will result in assigning the nearest p2 observation that occurs
    ## within N days. All p2 observations that occur within those N days are
    ## treated equally. If `max.dist` is non-null, p2 obs are only assigned to a
    ## p1 obs if they occur within `max.dist` (max.dist should be in same units
    ## as x,y coordinates specified in p1 and p2).
    ##
    ## With `unique = TRUE` each p2 obs is only assigned to one p1 obs,
    ## if multiple p1 obs are w/in the `max.dist` and `day` criteria then the
    ## p2 obs is assigned to the closest p1 obs. With `unique = FALSE`, an
    ## p2 obs is assigned to all p1 obs w/in the `max.dist` and `day` criteria,
    ## i.e., a single p2 obs can be assigned to multiple p1 obs.
    ##
    ## If both `max.days` and `max.dist` are null, the function assigns the
    ## closest p1 observation to each p1 observation, i.e., all p1 obs are
    ## assigned an ocean obs.
    ##
    ## p1 and p2 should be matrices with the first column giving the
    ## x-coordinate (e.g., longitude or easting) and the second column giving
    ## the y-coordinate (e.g., latitude or northing). If `max.days` is non-null,
    ## p1 and p2 should have a third column giving the julian day for the
    ## observation.
    ##
    ##
    ## p1 = 2 or 3 column matrix with columns lon,lat,jday
    ## p2 = 2 or 3 column matrix with columns lon,lat,jday
    ## method = method for calculating distance: "euclidean","lonlat"
    ## max.dist = max distance for matching
    ## max.days = max days for matching
    ## unique = logical, if TRUE only assign a p2 observation to the closest
    ##          p1 observation, if FALSE assign a p2 observation to all
    ##          p1 estimates w/in `max.days` and `max.dist`

    ## Safety checks
    if(!is.null(max.days)) {
        if(dim(p1)[2] < 3)
            stop("p1 matrix only has two columns")
        if(dim(p2)[2] < 3)
            stop("p2 matrix only has two columns")
    }

    p2.ind  <- 1:nrow(p2)
    ind.vec <- rep(NA, nrow(p1))
    dst.vec <- rep(NA, nrow(p1))
    if(!is.null(max.days)) p2_3 <- p2[ , 3]  ## move subset out of for loop

    
    per <- 0
    for(i in 1:nrow(p1)) {

        ## Account for `max.days`
        if(!is.null(max.days)) {
            day.range <- seq(p1[i , 3] - max.days, p1[i, 3] + max.days, 1)
            if(length(day.range) == 1) {
                p2.ind.day <- which(p2_3 == day.range)
            } else {
                p2.ind.day <- which(p2_3 >= min(day.range) &
                                    p2_3 <= max(day.range))
            }
            p2.day <- p2[p2.ind.day, ]
            if(length(p2.ind.day) == 1){
                p2.day <- matrix(p2.day, ncol = 3)
            }
        } else {
            p2.day <- p2
            p2.ind.day <- p2.ind
        }
        if(length(p2.ind.day) == 0) {
            next
        }

        ## Calculate distances
        if(method == "euclidean") {
            x.sq   <- (p1[i, 1] - p2.day[ , 1])^2
            y.sq   <- (p1[i, 2] - p2.day[ , 2])^2
            dist.i <- sqrt(x.sq + y.sq)
        }
        if(method == "lonlat") {
            ll1 <- p1[i, 1:2]
            ll2 <- p2.day[ , 1:2]
            dist.i <- geosphere::distHaversine(ll1, ll2)
        }
        min.ind <- which(dist.i == min(dist.i))[1]

        ## Account for `dist`
        if(is.null(max.dist) || dist.i[min.ind] <= max.dist) {
            ind.vec[i] <- p2.ind.day[min.ind]
            dst.vec[i] <- dist.i[min.ind]
        }
        
        if(loud && per != as.integer(i/nrow(p1) * 100)){
          per <- as.integer(i/nrow(p1) * 100)
          print(paste0(per, '% complete'))
        }

    }

    ## The assignment of p2 points to p1 points can result in the same p2
    ## point being assigned to multiple p1 points b/c the assignment
    ## is based only on days + distance.
    ##
    ## Make sure a single p2 point is only assigned to one p2 point.
    ## The criteria here is to only keep the assignment with the min distance
    ## b/w the p1 and p2 points.
    if(unique) {
        for(j in unique(ind.vec)) {

            ## get row indices for ind j
            rows.j <- which(!is.na(ind.vec) & ind.vec == j)

            ## get distances for ind j
            dist.j <- dst.vec[rows.j]

            ## get key w/ min distance
            row.min <- rows.j[which.min(dist.j)]

            ## get keys that should be set to NA
            rows.na <- rows.j[rows.j != row.min]

            ## set non-min-distance obs to NA
            ind.vec[rows.na] <- NA
            dst.vec[rows.na] <- NA
        }
    }

    out <- list(index = ind.vec,
                dist = dst.vec,
                max.dist = max.dist,
                max.days = max.days,
                unique = unique)
    return(out)
}

if(FALSE) {

    p1 <- matrix(1:20, ncol = 2)
    p2 <- matrix(11:30, ncol = 2)
    dist_min(p1, p2, method = "euclidean")
    dist_min(p1, p2, method = "euclidean", max.dist = 3)
    dist_min(p1, p2, method = "euclidean", max.dist = 10)
    dist_min(p1, p2, method = "euclidean", max.dist = 10, unique = TRUE)

    p1 <- matrix(c(1:20, rep(1, 10)), ncol = 3)
    p2 <- matrix(c(1:20, 1:10), ncol = 3)
    dist_min(p1, p2, method = "euclidean", max.days = 1)

    p1 <- matrix(c(1:20, rep(1, 10)), ncol = 3)
    p2 <- matrix(c(1:20, 3:12), ncol = 3)
    dist_min(p1, p2, method = "euclidean", max.days = 1)
    dist_min(p1, p2, method = "euclidean", max.days = 2)

}



## read_csv_dir --------------------------------------------
read_csv_dir <- function(path, pattern = "\\.csv$", verbose = TRUE) {
    ## Read all CSV in a directory
    ##
    ## Output is a data.frame of each file rbinded together
    ## non-existing columns in files are filled with NA
    ##
    ## path = path to directory with CSV files
    ## pattern = regex pattern to match for files
    ## verbose = print potentially useful information

    if(!dir.exists(path))
        stop("Path directory does not exist", call. = FALSE)

    fnames <- list.files(path, pattern = pattern, ignore.case = TRUE)
    out.lst <- vector("list", length(fnames))

    for(i in seq_along(fnames)) {
        i.name <- fnames[i]
        out.lst[[i]] <- data.table::fread(paste0(path, i.name), sep = ",",
                                          header = TRUE)
    }
    out.dt <- data.table::rbindlist(out.lst, fill = TRUE)
    out.df <- as.data.frame(out.dt)

    all.files   <- list.files(path, full.names = TRUE)
    n.all.files <- length(all.files[!file.info(all.files)$isdir])
    nfiles      <- length(fnames)

    if(verbose)
        cat("Read", nfiles, "of", n.all.files, "files", "\n")

    return(out.df)
}



## download_ne_bathymetry ----------------------------------
download_ne_bathymetry <- function(write.path, verbose = TRUE) {
  ## Download Natural Earth 10m bathymetric data
  ##
  ## write.path = path to dir to write out downloaded data
  
  if(verbose)
    cat("Downloading NE bathymetric data...", "\n")
  
  if(!dir.exists(write.path))
    dir.create(write.path, recursive = TRUE)
  
  url1 <- "https://www.naturalearthdata.com/"
  url2 <- "https://www.naturalearthdata.com/download/"
  url3 <- "10m/physical/ne_10m_bathymetry_all.zip"
  b.url <- paste0(url1, url2, url3)
  zip.file <- "./tmp_bath.zip"
  download.file(b.url, destfile = zip.file, method = "libcurl")
  unzip(zip.file, exdir = write.path)
  unlink(zip.file, recursive = TRUE)
}



## convert_coords ------------------------------------------
convert_coords <- function(data,
                           proj.orig,
                           proj.new,
                           xvar.orig = "lon",
                           yvar.orig = "lat",
                           xvar.new = "easting_m",
                           yvar.new = "northing_m",
                           plot = FALSE, ...) {
    ## Add transformed x,y coordinates to the input data.frame
    ##
    ## This function takes as input a data.frame with lon and lat columns
    ## specified by `xvar.orig` and `yvar.orig` and a projection string for the
    ## lon/lat (`proj.orig`) and (1) converts the coordinates to the projection
    ## given by `proj.new` and adds columns `xvar.new` and `yvar.new` to
    ## the input data.frame. The output of the function is the input `data` with
    ## two additional columns added whose names are specified by `xvar.new`
    ## and `yvar.new`.
    ##
    ## data = a data.frame with lon and lat columns
    ## proj.orig = string giving the projection of the lon/lat coords in the
    ##             input data
    ## proj.new = string giving the projection to convert to
    ## xvar.orig = longitude column name in `data`
    ## yvar.orig = latitude column name in `data`
    ## xvar.new = name of easting/lon converted column in output data.frame
    ## yvar.new = name of northing/lat converted column in output data.frame
    ## plot = logical, should the converted coords be plotted
    ## ... = passed to plot()

    if(xvar.new %in% names(data))
        stop("`xvar.new` already exists in data")
    if(yvar.new %in% names(data))
        stop("`yvar.new` already exists in data")

    if(!xvar.orig %in% names(data))
        stop("`xvar.orig` not found")
    if(!yvar.orig %in% names(data))
        stop("`yvar.orig` not found")

    dat <- data.table::copy(data)

    ## Add index column
    dat$index <- 1:nrow(dat)

    ## Subset lon and lat columns
    df.sub <- subset(dat, select = c(xvar.orig, yvar.orig, "index"))

    ## Convert to sf
    sf.sub <- sf::st_as_sf(df.sub, coords = c(1, 2), crs = proj.orig)

    ## Convert from proj.orig to proj.new
    sf.proj <- sf::st_transform(sf.sub, proj.new)
    ind     <- st_drop_geometry(sf.proj)
    coord   <- as.data.frame(st_coordinates(sf.proj))
    df.proj <- data.frame(index = ind$index,
                          x = coord$X,
                          y = coord$Y)
    names(df.proj) <- c("index", xvar.orig, yvar.orig)

    names(df.proj)[names(df.proj) == xvar.orig] <- xvar.new
    names(df.proj)[names(df.proj) == yvar.orig] <- yvar.new

    ## Add eastings and northings to input df
    df.out <- merge(x = dat, y = df.proj, by = "index")
    df.out <- df.out[order(df.out$index), ]
    df.out[["index"]] <- NULL

    if(plot) {
        plot(df.out[[xvar.new]], df.out[[yvar.new]], ...)
    }

    return(df.out)
}



## check_data ----------------------------------------------
check_data <- function(current,
                       target = NULL,
                       check.row.names = FALSE,
                       verbose = FALSE) {
    ## Perform sanity data checks
    ##
    ## The function returns TRUE if all tests pass and an error otherwise. If
    ## `target` is supplied the `all.equal` function is used to compare
    ## `current` and `target` using only the columns that `intersect`.
    ##
    ## current = data.frame to check
    ## target = optional data.frame to compare against
    ## check.row.names = logical, should row.names be checked
    ## verbose = logical, print potentially useful info

    if("read_order" %in% names(current)) {
        if(verbose) cat("Checking `read_order`...", "\n")

        tst <- sum(duplicated(current$read_order))
        if(tst != 0)
            stop("Duplicates in `read_order` column")

        tst1 <- current$read_order
        tst2 <- current$read_order[order(current$read_order)]
        if(any(tst1 != tst2))
            stop("`read_order` column not ordered")
    }

    if("temp_c" %in% names(current)) {
        if(verbose) cat("Checking `temp_c`...", "\n")

        if(any(current$temp_c < 0))
            stop("Temperature values below 0 found")
        if(any(current$temp_c > 30))
            stop("Temperature values above 30 found")
    }

    if("temp_100" %in% names(current)) {
        if(verbose) cat("Checking `temp_100`...", "\n")

        if(any(current$temp_c < 5))
            stop("Temperature values below 5 found")
        if(any(current$temp_c > 30))
            stop("Temperature values above 20 found")
    }

    if("wgt_total" %in% names(current)) {
        if(verbose) cat("Checking `wgt_total`...", "\n")

        if(any(current$wgt_total < 0))
            stop("Negative values in `wgt_total` column")
    }

    if("lat" %in% names(current)) {
        if(verbose) cat("Checking `lat`...", "\n")

        na <- sum(is.na(current$lat))
        if(na != 0)
            stop(na, "NAs in `lat` column")
    }

    if("date" %in% names(current)) {
        if(verbose) cat("Checking `date`...", "\n")

        na <- sum(is.na(current$date))
        if(na != 0)
            stop(na, "NAs in `date` column")

        yr <- as.numeric(format(current$date, format = "%Y"))
        if(any(current$year != yr))
            stop("`year` doesn't match year in `date`")

        if(!is(current$date, "Date"))
            stop("`Date` not of class 'Date'")
    }

    if("jday" %in% names(current)) {
        if(verbose) cat("Checking `jday`...", "\n")

        if(any(current$jday > 300))
            stop("jday values above 300 found")

        if(any(current$jday < 100))
            stop("jday values below 100 found")
    }

    if(check.row.names) {
        if(verbose) cat("Checking row.names...", "\n")
        c.rn <- as.numeric(row.names(current))

        d <- current
        row.names(d) <- NULL
        d.rn <- as.numeric(row.names(d))

        if(any(c.rn != d.rn))
            stop("row.names not set to NULL")
    }

    ae <- TRUE
    if(!is.null(target)) {
        if(verbose) cat("Comparing to target...", "\n")
        ae <- all_equal_intersect(target, current)
    }

    if(!isTRUE(ae)) {
        if(verbose) {
            print(ae)
        }
        stop("Error in hake data processing")
    }

    if(verbose) cat("Done!", "\n")
    return(ae)
}



## all_equal_intersect -------------------------------------
all_equal_intersect <- function(target, current, ...) {
    ## Call all.equal() using only intersecting columns
    ##
    ## target = a data.frame
    ## current = a data.frame
    ## ... = passed to all.equal()

    if(!is(target, "data.frame") || !is(current, "data.frame"))
        stop("Input is not class `data.frame`", call. = FALSE)

    c.by <- intersect(names(target), names(current))
    d1 <- subset(target, select = c.by)
    d2 <- subset(current, select = c.by)
    ae <- all.equal(d1, d2, ...)
    return(ae)
}

if(FALSE) {

    all_equal_intersect(iris, iris)
    all_equal_intersect(iris, iris[-1, ])
    all_equal_intersect(iris, as.matrix(iris))

}



## write_data_csv ------------------------------------------
write_data_csv <- function(data, write.path, verbose = TRUE) {
    ## Write input data to a CSV file
    ##
    ## This function writes `data` to a file with the same name as the `data`
    ## object (with dots replaced by underscores) as a CSV file in the directory
    ## specified by `write.path`. If the `write.path` directory doesn't exist,
    ## it is created
    ##
    ## data = data.frame with hake biomass data
    ## write.path = path to dir that CSV file is written to

    if(!dir.exists(write.path))
        dir.create(write.path)

    obj <- deparse(substitute(data))
    f.name <- gsub(".", "_", obj, fixed = TRUE)
    f.path <- paste0(write.path, f.name, ".csv")

    if(verbose)
        cat("Writing", f.path, "\n")

    write.csv(data, f.path, row.names = FALSE)
}



## <- --------------------------------------------
theme_simple <- function(base_size = 11,
                         base_family = "",
                         grid = FALSE,
                         mode = "light") {

    half_line <- base_size / 2
    if(mode == "dark") {
        grey0 <- "grey70"
        grey1 <- "grey30"
        grey2 <- "grey18"
        bg    <- "black"
    } else {
        grey0 <- "grey30"
        grey1 <- "grey70"
        grey2 <- "grey87"
        bg    <- NA
    }

    x <- theme_light(base_size = base_size, base_family = base_family) +
        theme(plot.background = element_rect(fill = bg, colour = NA),
              panel.background = element_rect(fill = bg, colour = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.ticks.length = unit(half_line / 2.2, "pt"),
              strip.background = element_rect(fill = NA, colour = NA),
              strip.text.x = element_text(colour = grey0),
              strip.text.y = element_text(colour = grey0),
              axis.text = element_text(colour = grey0),
              axis.title = element_text(colour = grey0),
              legend.title = element_text(colour = grey0, size = rel(0.9)),
              panel.border = element_rect(fill = NA, colour = grey1, size = rel(1.1)),
              legend.key.size = unit(0.9, "lines"),
              legend.text = element_text(size = rel(0.7), colour = grey0),
              legend.key = element_rect(colour = NA, fill = bg),
              legend.background = element_rect(colour = NA, fill = bg),
              plot.title = element_text(colour = grey0, size = rel(1)),
              plot.subtitle = element_text(colour = grey0, size = rel(.85)),
              plot.caption = element_text(colour = grey0, size = rel(0.8)),
              plot.tag = element_text(colour = grey0, size = rel(1.1)))
    if(grid) {
        x <- x + theme(panel.grid.major = element_line(color = grey2, size = 0.2))
    }
    return(x)
}

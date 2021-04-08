function (X, ..., r = NULL, rmax = NULL, breaks = NULL, correction = c("border", 
    "isotropic", "Ripley", "translate"), nlarge = 3000, 
    domain = NULL, var.approx = FALSE, ratio = FALSE) 
{
    verifyclass(X, "ppp")
    nlarge.given <- !missing(nlarge) && !is.null(nlarge)
    rfixed <- !is.null(r) || !is.null(breaks)
    npts <- npoints(X) #Done
    W <- X$window #Done
    areaW <- area(W) #Done
    lambda <- npts/areaW #Done
    lambda2 <- (npts * (npts - 1))/(areaW^2) #Done
    if (!is.null(domain)) {
        domain <- as.owin(domain)
        if (!is.subset.owin(domain, W)) 
            stop(paste(dQuote("domain"), "is not a subset of the window of X"))
        indom <- factor(inside.owin(X$x, X$y, domain), levels = c(FALSE, 
            TRUE))
        Kd <- Kdot(X %mark% indom, i = "TRUE", r = r, breaks = breaks, 
            correction = correction, ratio = ratio)
        Kd <- rebadge.fv(Kd, quote(K(r)), "K")
        return(Kd)
    }
    rmaxdefault <- rmax %orifnull% rmax.rule("K", W, lambda)
    if (is.infinite(rmaxdefault)) 
        rmaxdefault <- diameter(W)
    breaks <- handle.r.b.args(r, breaks, W, rmaxdefault = rmaxdefault)
    r <- breaks$r
    rmax <- breaks$max
    correction.given <- !missing(correction) && !is.null(correction)
    if (is.null(correction)) 
        correction <- c("border", "isotropic", "Ripley", 
            "translate")
    correction <- pickoption("correction", correction, 
        c(none = "none", border = "border", bord.modif = "bord.modif", 
            isotropic = "isotropic", Ripley = "isotropic", 
            trans = "translate", translate = "translate", 
            translation = "translate", rigid = "rigid", 
            good = "good", best = "best"), multi = TRUE)
    if ("good" %in% correction) 
        correction[correction == "good"] <- good.correction.K(X)
    correction <- implemented.for.K(correction, W$type, correction.given)
    alim <- c(0, min(rmax, rmaxdefault))
    can.do.fast <- breaks$even
    large.n <- (npts >= nlarge)
    large.n.trigger <- large.n && !correction.given
    fastcorrections <- c("border", "bord.modif", 
        "none")
    fastdefault <- "border"
    correction.fast <- all(correction %in% fastcorrections)
    will.do.fast <- can.do.fast && (correction.fast || large.n.trigger)
    asked <- correction.fast || (nlarge.given && large.n.trigger)
    if (asked && !can.do.fast) 
        warning("r values not evenly spaced - cannot use efficient code")
    if (will.do.fast) {
        ok <- correction %in% fastcorrections
        correction <- if (any(ok)) 
            correction[ok]
        else fastdefault
        bord <- any(correction %in% c("border", "bord.modif"))
        none <- any(correction == "none")
        if (!all(ok)) {
            corx <- c(if (bord) "border correction estimate" else NULL, 
                if (none) "uncorrected estimate" else NULL)
            corx <- paste(corx, collapse = " and ")
            message(paste("number of data points exceeds", 
                nlarge, "- computing", corx, "only"))
        }
        if (!rfixed) 
            r <- seq(from = 0, to = alim[2], length.out = length(r))
        if (bord) 
            Kb <- Kborder.engine(X, max(r), length(r), correction, 
                ratio = ratio)
        if (none) 
            Kn <- Knone.engine(X, max(r), length(r), ratio = ratio)
        if (bord && none) {
            Kn <- Kn[, names(Kn) != "theo"]
            yn <- fvnames(Kb, ".y")
            Kbn <- if (!ratio) 
                bind.fv(Kb, Kn, preferred = yn)
            else bind.ratfv(Kb, Kn, preferred = yn)
            return(Kbn)
        }
        if (bord) 
            return(Kb)
        if (none) 
            return(Kn)
    }
    do.fast.rectangle <- can.do.fast && is.rectangle(W) && spatstat.options("use.Krect") && 
        !any(correction == "rigid")
    if (do.fast.rectangle) {
        K <- Krect.engine(X, rmax, length(r), correction, ratio = ratio)
        attr(K, "alim") <- alim
    }
    else {
        Kdf <- data.frame(r = r, theo = pi * r^2)
        desc <- c("distance argument r", "theoretical Poisson %s")
        denom <- lambda2 * areaW
        K <- ratfv(Kdf, NULL, denom, "r", quote(K(r)), 
            "theo", NULL, alim, c("r", "%s[pois](r)"), 
            desc, fname = "K", ratio = ratio)
        rmax <- max(r)
        what <- if (any(correction %in% c("translate", 
            "isotropic"))) 
            "all"
        else "ijd"
        close <- closepairs(X, rmax, what = what)
        DIJ <- close$d
        gW <- NULL
        if (any(correction %in% c("translate", "rigid", 
            "isotropic"))) 
            gW <- setcov(W)
        if (any(correction == "none")) {
            wh <- whist(DIJ, breaks$val)
            numKun <- cumsum(wh)
            denKun <- lambda2 * areaW
            K <- bind.ratfv(K, data.frame(un = numKun), denKun, 
                "hat(%s)[un](r)", "uncorrected estimate of %s", 
                "un", ratio = ratio)
        }
        if (any(correction == "border" | correction == 
            "bord.modif")) {
            b <- bdist.points(X)
            I <- close$i
            bI <- b[I]
            RS <- Kount(DIJ, bI, b, breaks)
            if (any(correction == "bord.modif")) {
                denom.area <- eroded.areas(W, r)
                numKbm <- RS$numerator
                denKbm <- lambda2 * denom.area
                K <- bind.ratfv(K, data.frame(bord.modif = numKbm), 
                  data.frame(bord.modif = denKbm), "hat(%s)[bordm](r)", 
                  "modified border-corrected estimate of %s", 
                  "bord.modif", ratio = ratio)
            }
            if (any(correction == "border")) {
                numKb <- RS$numerator
                denKb <- lambda * RS$denom.count
                K <- bind.ratfv(K, data.frame(border = numKb), 
                  data.frame(border = denKb), "hat(%s)[bord](r)", 
                  "border-corrected estimate of %s", "border", 
                  ratio = ratio)
            }
        }
        if (any(correction == "translate")) {
            edgewt <- edge.Trans(dx = close$dx, dy = close$dy, 
                W = W, paired = TRUE, gW = gW, give.rmax = TRUE)
            wh <- whist(DIJ, breaks$val, edgewt)
            numKtrans <- cumsum(wh)
            denKtrans <- lambda2 * areaW
            h <- attr(edgewt, "rmax")
            numKtrans[r >= h] <- NA
            K <- bind.ratfv(K, data.frame(trans = numKtrans), 
                denKtrans, "hat(%s)[trans](r)", "translation-corrected estimate of %s", 
                "trans", ratio = ratio)
        }
        if (any(correction == "rigid")) {
            CW <- rotmean(gW)
            edgewt <- areaW/as.function(CW)(DIJ)
            wh <- whist(DIJ, breaks$val, edgewt)
            numKrigid <- cumsum(wh)
            denKrigid <- lambda2 * areaW
            h <- rmax.Rigid(X, gW)
            numKrigid[r >= h] <- NA
            K <- bind.ratfv(K, data.frame(rigid = numKrigid), 
                denKrigid, "hat(%s)[rigid](r)", "rigid motion-corrected estimate of %s", 
                "rigid", ratio = ratio)
        }
        if (any(correction == "isotropic")) {
            XI <- ppp(close$xi, close$yi, window = W, check = FALSE)
            edgewt <- edge.Ripley(XI, matrix(DIJ, ncol = 1))
            wh <- whist(DIJ, breaks$val, edgewt)
            numKiso <- cumsum(wh)
            denKiso <- lambda2 * areaW
            h <- boundingradius(W)
            numKiso[r >= h] <- NA
            K <- bind.ratfv(K, data.frame(iso = numKiso), denKiso, 
                "hat(%s)[iso](r)", "Ripley isotropic correction estimate of %s", 
                "iso", ratio = ratio)
        }
    }
    if (var.approx) {
        A <- areaW
        P <- perimeter(W)
        n <- npts
        rip <- 2 * ((A/(n - 1))^2) * (pi * r^2/A + 0.96 * P * 
            r^3/A^2 + 0.13 * (n/A) * P * r^5/A^2)
        if (!ratio) {
            K <- bind.fv(K, data.frame(rip = rip), "vR(r)", 
                "Ripley approximation to var(%s) under CSR", 
                "iso")
        }
        else {
            den <- (n - 1)^2
            ripnum <- den * rip
            ripden <- rep.int(den, length(rip))
            K <- bind.ratfv(K, data.frame(rip = ripnum), data.frame(rip = ripden), 
                "vR(r)", "Ripley approximation to var(%s) under CSR", 
                "iso")
        }
        if (W$type == "rectangle") {
            a1r <- (0.21 * P * r^3 + 1.3 * r^4)/A^2
            a2r <- (0.24 * P * r^5 + 2.62 * r^6)/A^3
            br <- (pi * r^2/A) * (1 - pi * r^2/A) + (1.0716 * 
                P * r^3 + 2.2375 * r^4)/A^2
            ls <- (A^2) * (2 * br - a1r + (n - 2) * a2r)/(n * 
                (n - 1))
            if (!ratio) {
                K <- bind.fv(K, data.frame(ls = ls), "vLS(r)", 
                  "Lotwick-Silverman approx to var(%s) under CSR", 
                  "iso")
            }
            else {
                den <- n * (n - 1)
                lsnum <- ls * den
                lsden <- rep.int(den, length(ls))
                K <- bind.ratfv(K, data.frame(ls = lsnum), data.frame(ls = lsden), 
                  "vLS(r)", "Lotwick-Silverman approx to var(%s) under CSR", 
                  "iso")
            }
        }
    }
    formula(K) <- . ~ r
    nama <- rev(colnames(K))
    fvnames(K, ".") <- setdiff(nama, c("r", "rip", 
        "ls"))
    unitname(K) <- unitname(X)
    if (ratio) 
        K <- conform.ratfv(K)
    return(K)
}
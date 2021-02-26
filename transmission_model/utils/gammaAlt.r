# These have been taken from the EnvStats library (https://cran.r-project.org/web/packages/EnvStats/index.html)
pgammaAlt = function (q, mean, cv = 1, lower.tail = TRUE, log.p = FALSE) 
{
    shape <- cv^-2
    scale <- mean/shape
    stats::pgamma(q = q, shape = shape, scale = scale, lower.tail = lower.tail, 
        log.p = log.p)
}

dgammaAlt = function (x, mean, cv = 1, log = FALSE)
{
    shape <- cv^-2
    scale <- mean/shape
    stats::dgamma(x = x, shape = shape, scale = scale, log = log)
}
rgammaAlt = function (n, mean, cv = 1) 
{
    shape <- cv^-2
    scale <- mean/shape
    stats::rgamma(n = n, shape = shape, scale = scale)
}

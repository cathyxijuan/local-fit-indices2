load("rmsea.ci.orig.mod.missing.below.RData")
load("rmsea.ci.orig.mod.missing.above.RData")
load("cfi.ci.orig.mod.missing.below.RData")
load("cfi.ci.orig.mod.missing.above.RData")
load("srmr.ci.orig.mod.missing.below.RData")
load("srmr.ci.orig.mod.missing.above.RData")

rmsea200below <- rmsea.ci.orig.mod.missing.below$n200

apply(rmsea200below, 2, function(x) which.max(abs(x-0.05)))
apply(rmsea200below, 2, function(x) which.min(abs(x-0.05)))

rmsea200above <- rmsea.ci.orig.mod.missing.above$n200

apply(rmsea200above, 2, function(x) which.max(abs(x-0.05)))
apply(rmsea200above, 2, function(x) which.min(abs(x-0.05)))

cfi200below <- cfi.ci.orig.mod.missing.below$n200[c(2,4),]

apply(cfi200below, 2, function(x) which.max(abs(x-0.05)))
apply(cfi200below, 2, function(x) which.min(abs(x-0.05)))

cfi200above <- cfi.ci.orig.mod.missing.above$n200[c(2,4),]

apply(cfi200above, 2, function(x) which.max(abs(x-0.05)))
apply(cfi200above, 2, function(x) which.min(abs(x-0.05)))



srmr200below <- srmr.ci.orig.mod.missing.below$n200[c(1, 2, 4,6),]

apply(srmr200below, 2, function(x) which.max(abs(x-0.05)))
apply(srmr200below, 2, function(x) which.min(abs(x-0.05)))

srmr200above <- srmr.ci.orig.mod.missing.above$n200[c(1, 2, 4,6),]

apply(srmr200above, 2, function(x) which.max(abs(x-0.05)))
apply(srmr200above, 2, function(x) which.min(abs(x-0.05)))





rmsea500below <- rmsea.ci.orig.mod.missing.below$n500

apply(rmsea500below, 2, function(x) which.max(abs(x-0.05)))
apply(rmsea500below, 2, function(x) which.min(abs(x-0.05)))

rmsea500above <- rmsea.ci.orig.mod.missing.above$n500

apply(rmsea500above, 2, function(x) which.max(abs(x-0.05)))
apply(rmsea500above, 2, function(x) which.min(abs(x-0.05)))

cfi500below <- cfi.ci.orig.mod.missing.below$n500[c(2,4),]

apply(cfi500below, 2, function(x) which.max(abs(x-0.05)))
apply(cfi500below, 2, function(x) which.min(abs(x-0.05)))

cfi500above <- cfi.ci.orig.mod.missing.above$n500[c(2,4),]

apply(cfi500above, 2, function(x) which.max(abs(x-0.05)))
apply(cfi500above, 2, function(x) which.min(abs(x-0.05)))



srmr500below <- srmr.ci.orig.mod.missing.below$n500[c(1, 2, 4,6),]


apply(srmr500below, 2, function(x) which.max(abs(x-0.05)))
apply(srmr500below, 2, function(x) which.min(abs(x-0.05)))

srmr500above <- srmr.ci.orig.mod.missing.above$n500[c(1, 2, 4,6),]


apply(srmr500above, 2, function(x) which.max(abs(x-0.05)))
apply(srmr500above, 2, function(x) which.min(abs(x-0.05)))


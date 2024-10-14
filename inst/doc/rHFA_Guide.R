## -----------------------------------------------------------------------------
library(rHFA)

in_corn = system.file('extdata', 'illinois_corn_yield.csv', package = "rHFA")
corn = read.csv(in_corn)

# convert columns to factors
facts = c('YEAR', 'COMPANY', 'HYBRID', 'SITE')
corn[facts] = lapply(corn[facts], as.factor)

head(corn, 10)
corn_full = corn

## -----------------------------------------------------------------------------
str(corn)

## -----------------------------------------------------------------------------
count_instances(corn, 'SITE', 'YEAR')

## -----------------------------------------------------------------------------
count_instances(corn, 'HYBRID', 'YEAR', TRUE) |>
  summary()

## -----------------------------------------------------------------------------
corn = filter_instances(
  data=corn,
  group_col='HYBRID',
  rep_col='YEAR',
  min_times=2,        # default is 2
  na.rm=TRUE          # ignore missing data instances
)

count_instances(corn, 'HYBRID', 'YEAR', TRUE) |>
  summary()

## -----------------------------------------------------------------------------
corn = filter_instances(
  data=corn,
  group_col='SITE',
  rep_col='YEAR',
  min_times=2,
  na.rm=TRUE
)

count_instances(corn, 'SITE', 'YEAR', na.rm=TRUE) |>
  summary()

## -----------------------------------------------------------------------------
corn$SITEYEAR = with(corn, 
                     paste(SITE, YEAR, sep='_'))

corn = filter_instances(corn,
                        group_col='SITEYEAR',
                        rep_col='HYBRID',
                        min_times=2,
                        na.rm=TRUE)

## -----------------------------------------------------------------------------
test = autofilter_instances(corn, 'SITE', 'YEAR', 'HYBRID', min_times=3, max_cycles=999)

## -----------------------------------------------------------------------------
pophfa = permute_hfa(corn,
                     level='population',  # default
                     site='SITE',
                     year='YEAR',
                     geno='HYBRID',
                     pheno='YIELD',
                     population=NA, # specify to run separately for each sub-population.
                     times=9, 
                     method='blup',
                     parallel=FALSE, 
                     seed=73956)

## -----------------------------------------------------------------------------
pophfa$home_field

## -----------------------------------------------------------------------------
pophfa$perms

## -----------------------------------------------------------------------------
pophfa$formula

## -----------------------------------------------------------------------------
head(pophfa$data)

## -----------------------------------------------------------------------------
mod = lm(pophfa$formula, 
         data=pophfa$data)
ss = get_ss(mod)
ss

## -----------------------------------------------------------------------------
get_home_site(pophfa$data, 
              geno='HYBRID', 
              site='SITE') |> 
  head()

## ----message=FALSE, warning=FALSE---------------------------------------------
corn_home = id_home(
  data=corn,
  site='SITE',
  year='YEAR', 
  geno='HYBRID', 
  pheno='YIELD',
  method='blup',     # need to update this to a switch to include median.
  verbose=FALSE
) 
subset(corn_home, HYBRID == corn[110, 'HYBRID']) |>
  head(10)

## -----------------------------------------------------------------------------
aggregate(rel_YIELD ~ SITE + HYBRID, 
          data=corn_home,
          subset=is_home,
          # data=subset(corn_home, is_home), 
          FUN=mean) |>
  summary()

## -----------------------------------------------------------------------------
genohfa = permute_hfa(corn,
                   level='genotype',
                   site='SITE',
                   year='YEAR',
                   geno='HYBRID',
                   pheno='YIELD',
                   times=9, 
                   method='blup',
                   parallel=FALSE, 
                   seed=73956)

## -----------------------------------------------------------------------------
genohfa_test = genohfa$home_field
head(genohfa_test)

## -----------------------------------------------------------------------------
tt = genohfa_test$p_val < 0.05
genohfa_test[tt, ]

## -----------------------------------------------------------------------------
genohfa_test$p_adj = p.adjust(genohfa_test$p_val, method='fdr')
tt = genohfa_test$p_val < 0.05
genohfa_test[tt, ]

## -----------------------------------------------------------------------------
geno_hfa = specialist_test(genohfa) |>
  generalist_test()
head(geno_hfa$home_field)

## -----------------------------------------------------------------------------
sitehfa = permute_hfa(corn,
                   level='site',
                   site='SITE',
                   year='YEAR',
                   geno='HYBRID',
                   pheno='YIELD',
                   times=9, 
                   method='blup',
                   parallel=FALSE, 
                   seed=73956)

## -----------------------------------------------------------------------------
sitehfa$home_field

## -----------------------------------------------------------------------------
yearhfa = permute_hfa(corn,
                   level='year',
                   site='SITE',
                   year='YEAR',
                   geno='HYBRID',
                   pheno='YIELD',
                   times=9, 
                   method='blup',
                   parallel=FALSE, 
                   seed=73956)

## -----------------------------------------------------------------------------
yearhfa$home_field

## -----------------------------------------------------------------------------
yhfa = yearhfa$home_field
plot(as.numeric(yhfa$YEAR), 
     yhfa$observed,
     xlab='Year',
     ylab='HFA')
# points(as.numeric(yhfa$YEAR), 
#        yhfa$observed)
points(as.numeric(yhfa$YEAR),
       yhfa$observed - yhfa$diff_q50, col='red')
legend('topright',
       c('expected', 'observed'),
       col=c('black', 'red'),
       pch=1)

## -----------------------------------------------------------------------------
lm(observed ~ as.numeric(YEAR), 
   data=yearhfa$home_field) |>
   summary()

## -----------------------------------------------------------------------------

corn_sites = system.file('extdata', 'illinois_sites.csv', package = "rHFA")
corn_sites = read.csv(corn_sites)

names(corn_sites) = gsub('City', 'SITE', names(corn_sites))
corn_sites[, 'SITE'] = tolower(corn_sites[, 'SITE'])

# home_dist = home_distance_spatial(
#   data=corn,
#   locations=corn_sites,
#   geno='HYBRID',
#   site='SITE',
#   lat='Lat',
#   long='Lon'
# )
# 
# home_dist[1:5, 1:5]

## -----------------------------------------------------------------------------
# envi_dist = home_distance_environmental(
#   data=corn,
#   locations=corn_sites,
#   geno='HYBRID',
#   site='SITE',
#   vars=c('bio_01', 'bio_17', 'BLDFIE_M_sl1_1km_ll', 'SLTPPT_M_sl1_1km_ll'),
#   scale=TRUE,
#   method='euclidean'
# )
# 
# envi_dist[1:5, 1:5]

## -----------------------------------------------------------------------------



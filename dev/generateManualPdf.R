library(subgroup.discovery)
path <- find.package("subgroup.discovery")
system(paste(shQuote(file.path(R.home("bin"), "R")), "CMD", "Rd2pdf", shQuote(path)))

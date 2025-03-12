# TODO: Fix: Do not hardcode the path...
files <- list.files("/Users/karcher/tactic/input/maps", full.names = TRUE)
files <- files[str_detect(files, "library_")]
files_modified <- map(files, \(x) {
    plate_number <- str_split(x, "_plate_")[[1]][2]
    plate_number <- str_replace(plate_number, ".csv", "")
    plate_number <- as.numeric(plate_number)
    if (
        plate_number %in% c(1, 2, 3)
    ) {
        xx <- read_csv(x) %>%
            mutate(biorep96 = 1)
    } else if (
        plate_number %in% c(4, 5, 6)
    ) {
        xx <- read_csv(x) %>%
            mutate(biorep96 = 2)
    } else {
        stop("This should never be reached - check your maps files for consistency")
    }

    xx <- xx %>%
        mutate(plt96 = plate_number)
    return(xx)
})
map2(files, files_modified, \(x, y) {
    write_csv(y, x)
})

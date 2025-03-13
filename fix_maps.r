# TODO: Fix: Do not hardcode the path...
files <- list.files("/Users/karcher/tactic/input/maps", full.names = TRUE)
files <- files[str_detect(files, "library_")]
files_modified <- map(files, \(x) {
    plate_number <- str_split(x, "_plate_")[[1]][2]
    plate_number <- str_replace(plate_number, ".csv", "")
    plate_number <- as.numeric(plate_number)

    if (plate_number %in% c(1, 2, 3)) {
        xx <- read_csv(x) %>%
            mutate(biorep96 = plate_number)
    } else if (plate_number == 4) {
        xx <- read_csv(x)
        xx$biorep96 <- NA
        xx$biorep96[1:22] <- 1
        xx$biorep96[23:44] <- 2
        xx$biorep96[45:66] <- 3
        xx$biorep96[67:88] <- 4
        xx$biorep96[89:96] <- 4 # These are some repeated genes...
    } else if (plate_number == 5) {
        xx <- read_csv(x) %>%
            mutate(biorep96 = 1)
    } else if (plate_number == 6) {
        xx <- read_csv(x) %>%
            mutate(biorep96 = 2)
    } else if (plate_number == 7) {
        xx <- read_csv(x) %>%
            mutate(biorep96 = 1)
    } else if (plate_number == 8) {
        xx <- read_csv(x) %>%
            mutate(biorep96 = 2)
    }
    xx <- xx %>%
        mutate(plt96 = plate_number)
    return(xx)
})
walk2(files, files_modified, \(x, y) {
    write_csv(y, x)
})

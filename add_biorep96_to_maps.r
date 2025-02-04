files <- list.files("/Users/karcher/tactic/input/maps", full.names = TRUE)
files <- files[str_detect(files, "library_")]
files_modified <- map(files, \(x) {
    # browser()
    if (
        str_detect(x, "_1") | str_detect(x, "_2") | str_detect(x, "_3")
    ) {
        read_csv(x) %>%
            mutate(biorep96 = 1)
    } else if (
        str_detect(x, "_4") | str_detect(x, "_5") | str_detect(x, "_6")) {
        read_csv(x) %>%
            mutate(biorep96 = 2)
    }
})
map2(files, files_modified, \(x, y) {
    write_csv(y, x)
})

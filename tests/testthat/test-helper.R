test_that("extract_time_space works", {
  wtemp <- data.frame(
    date = c("1979-04-02", "1979-04-02"),
    temp_0 = c(0, 0),
    temp_1 = c(1, 1),
    stringsAsFactors = FALSE)
  expect_length(extract_time_space(wtemp),
                2)
})

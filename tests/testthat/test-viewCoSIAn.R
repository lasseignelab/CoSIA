test_that("viewCoSIAn returns converted_id slot", {
    obj <- make_cosian()
    result <- viewCoSIAn(obj, "converted_id")
    expect_s3_class(result, "data.frame")
})

test_that("viewCoSIAn returns gex slot", {
    obj <- make_cosian()
    result <- viewCoSIAn(obj, "gex")
    expect_s3_class(result, "data.frame")
})

test_that("viewCoSIAn returns metric slot", {
    obj <- make_cosian()
    result <- viewCoSIAn(obj, "metric")
    expect_s3_class(result, "data.frame")
})

test_that("viewCoSIAn errors on invalid slot name", {
    obj <- make_cosian()
    expect_error(viewCoSIAn(obj, "invalid_slot"))
})

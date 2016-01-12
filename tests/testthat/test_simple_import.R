library(GenomicRanges)
library(soGGi)
test_ranges <- GRanges(seqnames = rep("chr1", 5),
                       ranges = IRanges(start = c(100, 200, 400, 800, 1000),
                                        width = c(50, 100, 250, 100, 50)),
                       strand = c("+", "-", "+", "+", "-"))

res_ranges1 <- GRanges(seqnames = rep("chr1", 5),
                       ranges = IRanges(start = c(124, 249, 524, 849, 1024),
                                        width = rep(1, 5)),
                       strand = c("+", "-", "+", "+", "-"))

res_ranges2 <- test_ranges
start(res_ranges2) <- start(res_ranges2) - width(res_ranges2)
end(res_ranges2) <- end(res_ranges2) + width(res_ranges2)/2

res_ranges3 <- GRanges(seqnames = rep("chr1", 5),
                       ranges = IRanges(start = c(50, 200, 350, 750, 1000),
                                       width = c(100, 150, 300, 150, 100)),
                       strand = c("+", "-", "+", "+", "-"))

context("Testing make_ranges function for use in imports")

expect_equal(test_ranges, make_ranges(testRanges = test_ranges, style = "region"))

expect_warning(make_ranges(testRanges = test_ranges, style = "centered", expand = NULL))
expect_equal(res_ranges1, make_ranges(testRanges = test_ranges, style = "centered", expand = 1))

expect_error(make_ranges(testRanges = test_ranges, style = "flanked", flank = NULL))
expect_equal(res_ranges2, make_ranges(testRanges = test_ranges, style = "flanked", flank = "100%"))
expect_equal(res_ranges3, make_ranges(testRanges = test_ranges, style = "flanked", 
                                      flankUp = 50, flankDown = 0))

expect_error(make_ranges(testRanges = test_ranges, style = "kittens"))
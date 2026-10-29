# Base-R analogue of withr::local_seed(): seed the global RNG for the calling
# test or fixture function and restore the previous stream, or its absence, when
# that frame exits. Seeded tests then leave no RNG state for later tests.
local_test_seed = function(seed, frame = parent.frame()) {
  global = globalenv()
  previous = if (exists('.Random.seed', envir = global, inherits = FALSE))
    get('.Random.seed', envir = global, inherits = FALSE)
  restore = function() {
    if (!is.null(previous)) assign('.Random.seed', previous, envir = global)
    else if (exists('.Random.seed', envir = global, inherits = FALSE))
      rm('.Random.seed', envir = global)
  }
  # Register the restore in the caller's frame, as withr's defer() does; later
  # seeds in the same frame are undone first.
  do.call(on.exit, list(as.call(list(restore)), add = TRUE, after = FALSE), envir = frame)
  set.seed(seed)
}

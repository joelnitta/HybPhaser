# Docker `--user` argument for the current host user

On Linux the container process must run as the host user so that it can
write to (and create sub-directories in) bind-mounted host directories.
Docker Desktop on macOS and Windows remaps bind-mount permissions, so
`--user` is neither needed nor always safe there.

## Usage

``` r
.docker_user()
```

## Value

`"<uid>:<gid>"`, or `NULL` when it should not be set.

# Run a Command in the rhybphaser Docker Container

Internal function to execute commands in the rhybphaser Docker container

## Usage

``` r
.run_docker(cmd, volumes = NULL, image = "joelnitta/rhybphaser:latest", ...)
```

## Arguments

- cmd:

  Character vector; command to run

- volumes:

  Named character vector; host paths to mount (names are container
  paths)

- image:

  Docker image name

- ...:

  Additional arguments passed to system2

## Value

Exit code from Docker command

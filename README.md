# boinc-apps

Container-based BOINC apps for projects run with
[boinc-stack](https://github.com/MBlenn/boinc-stack). Each app is a directory under `apps/`
that `boinc-app` (part of the boinc-stack server image) installs into a project. App images
that need more than the stack's base image are built here by GitHub Actions and pushed to GHCR.

| App | Purpose | Client requirements |
|---|---|---|
| [gnina](docs/gnina.md) | GPU molecular docking test app (redocking, cross-docking, error handling) with an RMSD validator | Linux, NVIDIA GPU, Docker or Podman with the NVIDIA Container Toolkit |
| [gnina_cpu](docs/gnina.md) | the same test workload on the CPU only, no GPU reserved | Linux, Docker or Podman |

## Installing an app into a project

On the server, as the project user:

```bash
podman exec <project>-daemons boinc-app install gnina --from https://github.com/MBlenn/boinc-apps
podman exec <project>-daemons boinc-app submit gnina 20
podman exec <project>-daemons boinc-app list
```

`--from https://github.com/MBlenn/boinc-apps@<tag or commit>` pins a specific state.
`install` copies the app definition into the project, so later changes here only take effect
after another `install`. To deploy changed app files, raise `APP_VERSION` in `app.env`.

`gnina` and `gnina_cpu` are opt-in (`APP_BETA=1`): volunteers only get them after enabling
"Run test applications" in their project preferences. Offer `gpucheck` (shipped with
boinc-stack, `boinc-app install gpucheck`) as the first step for new volunteers.

Volunteer host setup (BOINC client 8.2+, Podman/Docker, NVIDIA Container Toolkit) is described in
[boinc-stack's docs/volunteer-setup.md](https://github.com/MBlenn/boinc-stack/blob/main/docs/volunteer-setup.md).

## Layout

```
apps/<app>/          app definition (format: boinc-stack docs/docker-apps.md, "App format")
  app.env            name, version, plan class, validator/assimilator command lines
  plan_class.xml     plan class added to the project
  files/             Dockerfile, main.sh, job.toml (sent to the clients)
  templates/         create_work input/output templates
  server/            scripts copied to the project's bin/ (validator)
  jobs.py            job generator for "boinc-app submit"
  testdata/          input files used by jobs.py
app-images/<app>/    Containerfile of the app image (built in CI)
.github/workflows/   image builds and smoke tests
```

## Container images

The workflow `build-app-images` builds `ghcr.io/mblenn/boinc-app-<app>:<version>` on top of
`ghcr.io/mblenn/boinc-app-base:trixie` (built in boinc-stack). The clients pull these images
anonymously, so each GHCR package must be set to **public** after its first build
(Package settings → Change visibility).

App Dockerfiles in `apps/<app>/files/` must not contain `RUN` steps: they are built on the
volunteer's host by the BOINC client service, where rootless Podman cannot run them.

## License

MIT, see [LICENSE](LICENSE). Test data and third-party binaries keep their own licenses
(see the `NOTICE` files next to them).

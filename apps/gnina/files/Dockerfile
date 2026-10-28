# Built by docker_wrapper on the volunteer's host. The slot directory is mounted at /app.
# The app image (gnina binary + CUDA/cuDNN runtime libraries) is built in CI from
# app-images/gnina/ and pulled here; no RUN steps (see docs/docker-apps.md in boinc-stack).
FROM ghcr.io/mblenn/boinc-app-gnina:1.3.3
WORKDIR /app
CMD ["sh", "./main.sh"]

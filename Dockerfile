# Dockerfile for Hugging Face Spaces (Docker SDK).
# HF Spaces expects the app to listen on port 7860 by default.

FROM python:3.11-slim

# Create a non-root user for the runtime (HF Spaces requirement / good practice).
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR $HOME/app

# Install deps first so layer caching kicks in when only code changes.
COPY --chown=user:user requirements.txt ./
RUN pip install --no-cache-dir --user -r requirements.txt

# Copy the package + webapp source.
COPY --chown=user:user sgrna_array/ ./sgrna_array/
COPY --chown=user:user webapp/ ./webapp/
COPY --chown=user:user pyproject.toml ./

# Install the package itself so `from sgrna_array import ...` resolves cleanly.
RUN pip install --no-cache-dir --user .

EXPOSE 7860

# Production WSGI server. Two workers handle bursty traffic comfortably; HF Spaces
# free tier provides 16 GB RAM and 2 vCPUs.
CMD ["gunicorn", "--bind", "0.0.0.0:7860", "--workers", "2", "--timeout", "60", "webapp.app:app"]

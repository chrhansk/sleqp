set -e

PYTHON="python3"

cd bindings/python/dist

$PYTHON -m pip install --force-reinstall ./sleqp-*.tar.*

$PYTHON -c "import sleqp"

set -e

cd build

SLEQP_OCT_PACKAGE=(sleqp_octave*.tar.gz)

octave --eval "pkg install ${SLEQP_OCT_PACKAGE}"

octave --eval "pkg load sleqp; sleqp.info()"

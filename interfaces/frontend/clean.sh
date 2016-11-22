# @Author: jaume.bonet
# @Date:   2016-10-31 11:53:13
# @Last Modified by:   jaume.bonet
# @Last Modified time: 2016-10-31 11:53:25

export VENDOR='app/vendor'

rm -rf $VENDOR/jquery
rm -rf $VENDOR/fontawesome/less $VENDOR/fontawesome/scss
rm -rf $VENDOR/bootstrap/fonts $VENDOR/bootstrap/grunt $VENDOR/bootstrap/js $VENDOR/bootstrap/less
rm -rf $VENDOR/bootstrap/dist/js

rm -f $VENDOR/*/.*
rm -f $VENDOR/*/*.json
rm -f $VENDOR/*/*.md
rm -f $VENDOR/*/LICENSE
rm -f $VENDOR/*/Gruntfile.js
rm -f $VENDOR/*/bower.* $VENDOR/*/package.*
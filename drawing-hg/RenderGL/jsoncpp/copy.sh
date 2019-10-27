rm -f libs/libjson.a
libfile=$(find libs -name 'lib*.a')
echo "cp $libfile libs/libjson.a"
cp -fp $libfile libs/libjson.a
# Generate the case
if [ ! -f 2H4E2H.absolute.yml ]; then
  # Create an Architecture defined case
  topo.case -name 2H4E2H -architecture 2H.4E.2H
  # Apply corrections to make an absolute case
  topo.absolute -case 2H4E2H.yml -corrections corrections.yml -caseout 2H4E2H.absolute
fi

protocol=protocols/${1}_protocol.yml

echo 'Requested protocol: ' ${protocol}

if [ -f ${protocol} ]; then
  topo.protocol -case 2H4E2H.absolute.yml -protocol ${protocol}
else
  echo 'Unkown protocol prefix.'
fi


while [[ $# > 0 ]]
do
  echo -n "$1: "
  zcat  "$1" | awk '
    BEGIN { x=0; on=0; }
    /^_struct_sheet_order.sense/ { on=1; next; }
    /^loop_|^#/ { on=0; next; }
    { if (on==1) {x=x+1;}  }
    END { print "Num of struct sheet order lines:", x;  }
  '
  shift
done

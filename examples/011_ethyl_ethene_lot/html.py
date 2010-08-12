

header = """<?xml version='1.0' encoding='UTF-8'?>
<html><head><title>%s</title>
<style type='text/css'>
body { font-family: sans; }
td { text-align: right; width: 50px; padding: 0px 10px; }
table, td, th, tr { border: solid 1px black; }
table { table-layout: fixed; border-collapse: collapse; page-break-after: always }
</style>
</head><body>"""

footer = "</body></html>"


def print_table(f, rows):
    print >> f, "<table style='border-color:black'>"
    for row in rows:
        print >> f, "<tr>%s</tr>" % ("".join(row).encode('UTF-8'))
    print >> f, "</table>"

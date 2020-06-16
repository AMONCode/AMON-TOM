from sshtunnel import SSHTunnelForwarder
import MySQLdb as db
import pandas as pd

# ssh variables
host = '3.13.26.235'
localhost = '127.0.0.1'
ssh_username = 'ubuntu'
ssh_private_key = '/var/www/amonTOM/.ssh/pair_of_keys.pem'

with open('/var/www/amonTOM/amonTOM/amon_db_password.txt') as f:
    AMON_DB_PASSWORD = f.read().strip()

# database variables
user='amon'
password=AMON_DB_PASSWORD
database='AMON_test2'


def query(q):
     with SSHTunnelForwarder(
          (host, 22),
          ssh_username=ssh_username,
          ssh_private_key=ssh_private_key,
          remote_bind_address=(localhost, 3306)
     ) as server:
          conn = db.connect(host=localhost,
          port=server.local_bind_port,
          user=user,
          passwd=password,
          db=database)
          return pd.read_sql_query(q, conn)

df = query('select * from alert limit 5;')
df.head()

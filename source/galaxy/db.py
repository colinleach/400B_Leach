import psycopg2
import yaml
from pathlib import Path
home = str(Path.home())

class DB():
    """
    A simple wrapper class for connecting to the PostgreSQL database.

    Takes no argmuments. Relies on having connection information in
    `~/dbconn.yaml`. 
    """

    def __init__(self):
        "Reads the connection parameters, makes the connection and a cursor"

        params = self.read_params()

        inf = f"dbname={params['dbname']} user={params['username']}"
        inf += f"  host='{params['host']}' password={params['password']}"
        self.connection = psycopg2.connect(inf)
        self.connection.autocommit = True
        self.cursor  = self.connection.cursor()

    def read_params(self):
        "Needs the yaml parameter file to be in the user's home directory"

        filename = Path.home() / 'dbconn.yaml'
        with open(filename) as file:
            params = yaml.full_load(file)
            return params

    def get_cursor(self):
        "A simple getter method"

        return self.cursor

    def run_query(self, query):
        """
        Runs a SQL query (typically SELECT)

        Returns results in Python list format (not numpy)
        """

        self.cursor.execute(query)
        return self.cursor.fetchall()
 

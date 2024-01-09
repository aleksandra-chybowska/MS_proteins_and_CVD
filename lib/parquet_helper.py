import pyarrow as pa
import pyarrow.parquet as pq


def write_parquet(df, path):
    table = pa.Table.from_pandas(df)
    pq.write_table(table, path)


def read_parquet(path):
    table = pq.read_table(path)
    return table.to_pandas()


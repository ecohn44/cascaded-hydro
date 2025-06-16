def deduplicate_data(df, int_freq):

    # Remove data not on the int_freq min interval
    if int_freq == 15:
        df = df[df.index.minute.isin([0, 15, 30, 45]) & (df.index.second == 0) & (df.index.microsecond == 0)]
    if int_freq == 30:
        df = df[df.index.minute.isin([0, 30]) & (df.index.second == 0) & (df.index.microsecond == 0)]

    old_len = len(df)
    
    # Identify all duplicated rows based on the datetime index
    duplicate_rows = df[df.index.duplicated(keep=False)]  # keep=False keeps all copies
    
    # Extract unique dates (without time) from duplicated timestamps
    duplicate_dates = np.unique(duplicate_rows.index.date)
    
    # Print results
    print("Duplicate Dates:")
    for date in duplicate_dates:
        print(f"  - {date}")
    print("\nDuplicated Data:")
    print(duplicate_rows)

    # Remove duplicates
    df = df[~df.index.duplicated(keep='first')]

    new_len = len(df)
    removed_rows = old_len - new_len
    
    print(f"Removed {removed_rows} rows")
    print()
    
    return df


def drop_nan(df):
    before = len(df)
    df_cleaned = df.dropna(subset=['00065'])
    after = len(df_cleaned)
    removed = before - after
    print(f"Removed {removed} rows with NaN in '00065'")
    print()
    return df_cleaned
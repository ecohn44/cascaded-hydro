import numpy as np
import pandas as pd
import seaborn as sns
from datetime import datetime
import matplotlib.pyplot as plt
import dataretrieval.nwis as nwis
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score



def deduplicate_data(df, int_freq):

    print_block = False

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
    if print_block: 
        print("Duplicate Dates:")
        for date in duplicate_dates:
            print(f"  - {date}")
        print("\nDuplicated Data:")
        print(duplicate_rows)

    # Remove duplicates
    df = df[~df.index.duplicated(keep='first')]

    new_len = len(df)
    removed_rows = old_len - new_len
    
    print(f"Removed {removed_rows} duplicate rows")
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

def impute_missing_data(df, int_freq):

    print_block = False

    freq=str(int_freq) + 'min'
    
    # Create full expected range at specified frequency
    full_range = pd.date_range(start=df.index.min(), end=df.index.max(), freq=freq)

    # Find missing timestamps
    missing = full_range.difference(df.index)

    # Print missing days report
    missing_days = pd.Series(missing.date).unique()
    
    if print_block:
        if len(missing_days) > 0:
            print("Missing data on the following days:")
            for day in missing_days:
                print(f"- {day}")
        else:
            print("No missing days.")

    imputed_count = 0

    # Impute missing values using average of the endpoints
    for missing_time in missing:
        # Find the two nearest data points surrounding the missing timestamp
        prev_point = df[df.index <= missing_time].iloc[-1:]  # Last point before the missing timestamp
        next_point = df[df.index > missing_time].iloc[:1]    # First point after the missing timestamp

        if not prev_point.empty and not next_point.empty:
            avg_value = (prev_point['00065'].values[0] + next_point['00065'].values[0]) / 2
            df.loc[missing_time] = avg_value
            imputed_count += 1
            if print_block:
                print(f"Imputed missing data at {missing_time} with value {avg_value:.2f}")
        elif not prev_point.empty:
            df.loc[missing_time] = prev_point['00065'].values[0]
            imputed_count += 1
            if print_block:
                print(f"Imputed missing data at {missing_time} with value {prev_point['00065'].values[0]:.2f}")
        elif not next_point.empty:
            df.loc[missing_time] = next_point['00065'].values[0]
            imputed_count += 1
            if print_block:
                print(f"Imputed missing data at {missing_time} with value {next_point['00065'].values[0]:.2f}")
        else:
            print(f"Unable to impute missing data at {missing_time}.")

    print(f"Imputed {imputed_count} rows")

    # Resort imputed data
    df = df.sort_index()

    return df

def create_date_range(start_date, end_date):
    start_date_obj = datetime.strptime(start_date, '%Y-%m-%d') 
    start_new_format = start_date_obj.strftime('%m-%d-%Y') 
    end_date_obj = datetime.strptime(end_date, '%Y-%m-%d') 
    end_new_format = end_date_obj.strftime('%m-%d-%Y') 
    date_range = start_new_format + ' to ' + end_new_format

    return date_range

def preprocess(site_no, int_freq, start_date, end_date):
    
    print(nwis.get_record(sites=site_no, service='site').station_nm[0])

    df = nwis.get_record(sites=site_no, service='iv', start=start_date, end=end_date)
    df.index = df.index.tz_convert('America/Los_Angeles')  # Convert to Pacific Time
    df.index = df.index.tz_localize(None) # Remove timezone info 
    df = df[['site_no','00065']]

    # Deduplicate data
    df = deduplicate_data(df, int_freq)

    # Remove NaNs
    df = drop_nan(df)

    # Impute missing data
    df = impute_missing_data(df, int_freq)

    return df

def dataset_merge(df1, df2, train_start_year, season):
    
    if season == "full":
        month_start = '01-01'
        month_end = '12-31'
    else:
        month_start = '04-01'
        month_end = '08-31'

    # Apply Seasonality
    train_start_date = train_start_year + '-' + month_start
    df_range = df1[(df1.index >= train_start_date)]

    mask = (
        (df_range.index.strftime('%m-%d') >= month_start) &
        (df_range.index.strftime('%m-%d') <= month_end)
    )

    if season == "dry":
        # isolate summer months
        df_range = df_range[mask]
    elif season == "wet":
        # isolate wet months
        df_range = df_range[~mask]
    else:
        # full year so no mask needed
        pass

    df_range = df_range.rename(columns={'00065': 'y'}) ## df: inflow downstream

    # Merge in upstream data
    df_range = df_range.merge(df2, left_index=True, right_index=True, how='inner')
    df_range = df_range.rename(columns={'00065': 'x'})
    df_range = df_range[['y', 'x']]

    # Encode Day of Year
    df_range['DoY'] = pd.to_datetime(df_range.index).strftime('%j').astype(int)
    df_range['day_sin'] = np.sin(2 * np.pi * df_range['DoY'] / 365)
    df_range['day_cos'] = np.cos(2 * np.pi * df_range['DoY'] / 365)

    # Normalize features
    scaler = StandardScaler() # unit variance & zero mean
    df_range_norm = scaler.fit_transform(df_range)
    df_range['y_norm'] = df_range_norm[:,0]
    df_range['x_norm'] = df_range_norm[:,1]

    return df_range

def create_lag_features(df_reg, p, up_feat, down_feat):
    # Adding upstream lag features to the DataFrame
    if up_feat: 
        for i in range(1, p):  # Creating lag features up to p time units
            df_reg[f'X_Lag_{i}'] = df_reg['x_norm'].shift(i)

    if down_feat: 
        # Adding downstream lag features to the DataFrame
        for i in range(1, p):  # Creating lag features up to p time units
            df_reg[f'Y_Lag_{i}'] = df_reg['y_norm'].shift(i)

    # Drop rows with NaN values resulting from creating lag features
    df_reg.dropna(inplace=True)
    df_reg = df_reg.rename(columns={1: '1'})

    # Define features and target
    if up_feat and down_feat: 
        feature_cols = ['day_sin', 'day_cos', 'x_norm'] + [f'X_Lag_{lag}' for lag in range(1, p)] + [f'Y_Lag_{lag}' for lag in range(1, p)]
    if up_feat and not down_feat: 
        feature_cols = ['day_sin', 'day_cos', 'x_norm'] + [f'X_Lag_{lag}' for lag in range(1, p)] 
    if down_feat and not up_feat: 
        feature_cols = ['day_sin', 'day_cos'] + [f'Y_Lag_{lag}' for lag in range(1, p)]
    
    return df_reg, feature_cols

def print_test_stats(y_test, y_pred):
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    mae = mean_absolute_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)

    print(f'RMSE: {rmse:.3f}')
    print(f'MAE: {mae:.3f}')
    print(f'R²: {r2:.3f}')

def plot_forecasts(model, y_test, y_pred):
    plt.figure(figsize=(10, 4))
    plt.plot(y_test.index, y_test.values, label='Actual')
    plt.plot(y_test.index, y_pred.values, label='Predicted', linestyle='--')
    plt.title(model + ': Model Predictions')
    plt.legend()
    plt.xlabel('Time Index')
    plt.ylabel('Gage Height [normalized]')
    plt.tight_layout()
    plt.show()

def plot_variance(model, fitted_vals, residuals, bounds=True):
    plt.figure(figsize=(10, 2.5))
    sns.scatterplot(x=fitted_vals, y=residuals, alpha=0.7)
    plt.axhline(0, color='red', linestyle='--')
    if bounds:
        plt.ylim([-.6, .6])
    plt.xlabel("Fitted values")
    plt.ylabel("Residuals")
    plt.title(model + ": Fitted Values vs. Residuals")
    plt.grid(True)
    plt.show()

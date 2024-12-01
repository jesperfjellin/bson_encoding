from pathlib import Path
import gzip
from bson import encode, decode
from bson.codec_options import CodecOptions
import json

def format_size(bytes_size):
    """
    Convert a file size in bytes to a human-readable string.
    
    Args:
        bytes_size (int): File size in bytes.
        
    Returns:
        str: Human-readable file size.
    """
    for unit in ['B', 'KB', 'MB', 'GB']:
        if bytes_size < 1024:
            return f"{bytes_size:.2f} {unit}"
        bytes_size /= 1024
    return f"{bytes_size:.2f} TB"

def delta_encode(coords, scale_factor=1e6):
    """
    Delta encode a list of coordinates after scaling them to integers.
    
    Args:
        coords (list): List of [x, y] or [x, y, z] coordinates.
        scale_factor (float): Factor to scale coordinates for integer conversion.
        
    Returns:
        list: Delta encoded coordinates.
    """
    # Scale coordinates and convert to integers
    scaled_coords = [[int(round(coord * scale_factor)) for coord in point] for point in coords]
    # Start by storing the first point in full
    encoded_coords = [scaled_coords[0]]
    for i in range(1, len(scaled_coords)):
        # Delta encode each subsequent coordinate
        delta = [scaled_coords[i][j] - scaled_coords[i-1][j] for j in range(len(scaled_coords[i]))]
        encoded_coords.append(delta)
    return encoded_coords

def delta_decode(encoded_coords, scale_factor=1e6, decimal_places=6):
    """
    Delta decode a list of coordinates and scale them back to floating points.
    
    Args:
        encoded_coords (list): Delta encoded coordinates.
        scale_factor (float): Factor used during encoding to scale coordinates.
        decimal_places (int): Number of decimal places to round the decoded coordinates.
        
    Returns:
        list: Decoded coordinates.
    """
    # Restore the original scaled coordinates from the encoded deltas
    decoded_scaled_coords = [encoded_coords[0]]
    for delta in encoded_coords[1:]:
        next_point = [decoded_scaled_coords[-1][j] + delta[j] for j in range(len(delta))]
        decoded_scaled_coords.append(next_point)
    # Convert back to floating-point coordinates and round
    coords = [[round(coord / scale_factor, decimal_places) for coord in point] for point in decoded_scaled_coords]
    return coords

def encode_geojson(file_path, scale_factor=1e6):
    """
    Encode a GeoJSON file using delta encoding and BSON compression.
    
    Args:
        file_path (str): Path to the input GeoJSON file.
        scale_factor (float): Factor to scale coordinates for integer conversion.
    """
    # Read the original GeoJSON file
    with open(file_path, 'r') as f:
        geojson_data = json.load(f)

    # Prepare a list to hold encoded features
    encoded_features = []

    # Loop through each feature in the GeoJSON
    for feature in geojson_data['features']:
        geom = feature['geometry']
        properties = feature.get('properties', {})
        geometry_type = geom['type']
        original_coords = geom['coordinates']

        if geometry_type == 'Point':
            # Scale and convert to integers
            scaled_coords = [int(round(coord * scale_factor)) for coord in original_coords]
            geom['coordinates'] = scaled_coords
            properties['encoding'] = 'scaled'
        elif geometry_type in ['LineString', 'MultiPoint']:
            # Delta encode the coordinates
            encoded_coords = delta_encode(original_coords, scale_factor)
            geom['coordinates'] = encoded_coords
            properties['encoding'] = 'delta'
        elif geometry_type == 'Polygon':
            # Delta encode each ring separately
            encoded_coords = []
            start_points = []
            for ring in original_coords:
                encoded_ring = delta_encode(ring, scale_factor)
                encoded_coords.append(encoded_ring)
                # Store the scaled first point of each ring
                scaled_ring0 = [int(round(coord * scale_factor)) for coord in ring[0]]
                start_points.append(scaled_ring0)
            geom['coordinates'] = encoded_coords
            properties['encoding'] = 'delta'
            properties['start_points'] = start_points  # Store start points for decoding
        elif geometry_type == 'MultiLineString':
            # Delta encode each LineString
            encoded_coords = []
            for line in original_coords:
                encoded_line = delta_encode(line, scale_factor)
                encoded_coords.append(encoded_line)
            geom['coordinates'] = encoded_coords
            properties['encoding'] = 'delta'
        elif geometry_type == 'MultiPolygon':
            # Delta encode each Polygon
            encoded_coords = []
            for polygon in original_coords:
                encoded_polygon = []
                for ring in polygon:
                    encoded_ring = delta_encode(ring, scale_factor)
                    encoded_polygon.append(encoded_ring)
                encoded_coords.append(encoded_polygon)
            geom['coordinates'] = encoded_coords
            properties['encoding'] = 'delta'
        else:
            # Unsupported geometry type
            continue

        # Update the feature
        feature['geometry'] = geom
        feature['properties'] = properties
        encoded_features.append(feature)

    # Prepare the encoded data
    encoded_data = {
        'type': 'FeatureCollection',
        'features': encoded_features
    }

    # Use Pathlib to handle file paths
    path = Path(file_path)
    encoded_file_path = path.parent / f"encoded_{path.stem}.bson.gz"

    # Serialize the data to BSON
    bson_data = encode(encoded_data)

    # Compress the BSON data with GZIP
    with gzip.open(encoded_file_path, 'wb') as f:
        f.write(bson_data)

    print(f"Encoded and compressed data saved as '{encoded_file_path}'")

def decode_geojson(file_path, scale_factor=1e6, decimal_places=6):
    """
    Decode a BSON-compressed GeoJSON file back to standard GeoJSON.
    
    Args:
        file_path (str): Path to the encoded BSON GZIP file.
        scale_factor (float): Factor used during encoding to scale coordinates.
        decimal_places (int): Number of decimal places to round the decoded coordinates.
    """
    # Read the encoded and compressed BSON file
    with gzip.open(file_path, 'rb') as f:
        bson_data = f.read()

    # Deserialize the BSON data
    encoded_data = decode(bson_data)

    # Prepare a list to hold decoded features
    decoded_features = []

    # Loop through each feature
    for feature in encoded_data['features']:
        geom = feature['geometry']
        properties = feature.get('properties', {})
        geometry_type = geom['type']
        encoded_coords = geom['coordinates']

        encoding = properties.get('encoding')

        if geometry_type == 'Point' and encoding == 'scaled':
            # Convert back to floating-point coordinates
            coords = [coord / scale_factor for coord in encoded_coords]
            geom['coordinates'] = coords
        elif encoding == 'delta':
            if geometry_type in ['LineString', 'MultiPoint']:
                # Delta decode the coordinates
                decoded_coords = delta_decode(encoded_coords, scale_factor, decimal_places)
                geom['coordinates'] = decoded_coords
            elif geometry_type == 'Polygon':
                # Delta decode each ring separately
                # Retrieve the list of start_points
                start_points = properties.get('start_points', [])
                if len(start_points) != len(encoded_coords):
                    print(f"Warning: Number of start points does not match number of rings for Polygon.")
                decoded_coords = []
                for ring, start_pt in zip(encoded_coords, start_points):
                    # start_pt is already scaled, no need to scale here
                    decoded_ring = delta_decode(ring, scale_factor, decimal_places)
                    decoded_coords.append(decoded_ring)
                geom['coordinates'] = decoded_coords
                # Remove start_points after decoding
                if 'start_points' in properties:
                    del properties['start_points']
            elif geometry_type == 'MultiLineString':
                # Delta decode each LineString
                decoded_coords = []
                for line in encoded_coords:
                    decoded_line = delta_decode(line, scale_factor, decimal_places)
                    decoded_coords.append(decoded_line)
                geom['coordinates'] = decoded_coords
            elif geometry_type == 'MultiPolygon':
                # Delta decode each Polygon
                decoded_coords = []
                for polygon in encoded_coords:
                    decoded_polygon = []
                    for ring in polygon:
                        decoded_ring = delta_decode(ring, scale_factor, decimal_places)
                        decoded_polygon.append(decoded_ring)
                    decoded_coords.append(decoded_polygon)
                geom['coordinates'] = decoded_coords
            else:
                # Unsupported geometry type
                pass
        else:
            # If not encoded, assume coordinates are already in floating-point form
            pass

        # Remove the 'encoding' property and any other encoding-related properties
        if 'encoding' in properties:
            del properties['encoding']
        if 'start_points' in properties:
            del properties['start_points']

        # Update the feature
        feature['geometry'] = geom
        feature['properties'] = properties
        decoded_features.append(feature)

    # Prepare the decoded data
    decoded_data = {
        'type': 'FeatureCollection',
        'features': decoded_features
    }

    # Use Pathlib to handle file paths
    path = Path(file_path)
    decoded_file_path = path.parent / f"decoded_{path.stem.replace('encoded_', '')}.geojson"

    # Write the decoded data to a GeoJSON file with compact separators
    with open(decoded_file_path, 'w') as f:
        json.dump(decoded_data, f, separators=(',', ':'))

    print(f"Decoded GeoJSON saved as '{decoded_file_path}'")

def print_file_sizes(original_file, encoded_file, decoded_file):
    """
    Print the sizes of the original, encoded, and decoded files.
    
    Args:
        original_file (Path): Path object for the original GeoJSON file.
        encoded_file (Path): Path object for the encoded BSON GZIP file.
        decoded_file (Path): Path object for the decoded GeoJSON file.
    """
    original_size = original_file.stat().st_size
    encoded_size = encoded_file.stat().st_size
    decoded_size = decoded_file.stat().st_size

    print("\nFile Sizes:")
    print(f"Original GeoJSON: {format_size(original_size)}")
    print(f"Encoded BSON GZIP: {format_size(encoded_size)}")
    print(f"Decoded GeoJSON: {format_size(decoded_size)}\n")

# Example usage:
if __name__ == "__main__":
    original_file_path = r'C:\Rust\Projects\rust_bindings\test_points.geojson'
    encoded_file_path = Path(original_file_path).parent / f"encoded_{Path(original_file_path).stem}.bson.gz"
    decoded_file_path = Path(original_file_path).parent / f"decoded_{Path(encoded_file_path).stem.replace('encoded_', '')}.geojson"

    # Encode the GeoJSON file
    encode_geojson(original_file_path)

    # Decode the encoded file
    decode_geojson(encoded_file_path)

    # Print the file sizes
    print_file_sizes(Path(original_file_path), encoded_file_path, decoded_file_path)

export type ProjectedPoint = [number, number];
export type Projection = (latitude: number, longitude: number) => ProjectedPoint | null;

const RAD = Math.PI / 180;

export function longitude180(longitude: number): number {
  return ((longitude + 180) % 360 + 360) % 360 - 180;
}

export function orthographic(
  width: number, height: number, latitude: number, longitude: number,
): { radius: number; project: Projection } {
  const radius = Math.min(width, height) * 0.44;
  const center = latitude * RAD;
  const sinCenter = Math.sin(center);
  const cosCenter = Math.cos(center);
  return {
    radius,
    project(lat, lon) {
      const angle = lat * RAD;
      const delta = longitude180(lon - longitude) * RAD;
      const visibility = sinCenter * Math.sin(angle) + cosCenter * Math.cos(angle) * Math.cos(delta);
      if (visibility < -1e-12) return null;
      return [
        width / 2 + radius * Math.cos(angle) * Math.sin(delta),
        height / 2 - radius * (cosCenter * Math.sin(angle) - sinCenter * Math.cos(angle) * Math.cos(delta)),
      ];
    },
  };
}

export function equirectangular(width: number, height: number): Projection {
  return (lat, lon) => [((lon + 180) / 360) * width, ((90 - lat) / 180) * height];
}
